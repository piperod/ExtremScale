#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

double mvallgather(const double *x, const double *mat, double *y,
		   int gn, int ln, int sr, int er, double *tcomm);
double mvring(const double *x, const double *mat, double *y,
	      int gn, int ln, int sr, int er, double *tcomm);
void mvmult(const double  x, const double  mat,
	    double t y, int rows, int cols, int colsdim);
void initVars(double *x, double *y, double *mat, int ln, int gn, int sr);
void printVec(double *y, int ln);

static int verbose = 0;

int main(int argc, char *argv[])
{
    int gn, ln, sr, er;
    int wsize, wrank, k;
    double *x, *y, *mat;
    double t_ag, t_ring, tcomm_ag, tcomm_ring, t, t2;
    const int ntest = 10;

    MPI_Init(0,0);

    /* Define the problem */
    gn  = 128;
    if (argc > 1 && argv[1]) gn = atoi(argv[1]);

    /* Define the decomposition */
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    ln  = gn / wsize;
    sr  = wrank * ln;
    er  = sr + ln - 1;
    if (ln * wsize != gn) {
	fprintf(stderr, "gn must be a multiple of wsize\n");
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Create the vectors and matrix */
    x = (double *)malloc(ln*sizeof(double));
    y = (double *)malloc(ln*sizeof(double));
    mat = (double *)malloc(ln*gn*sizeof(double));

    /* Perform the matrix-vector multiplication */
    t_ag = 1e6; t_ring = 1e6;
    tcomm_ag = 1e6; tcomm_ring = 1e6;
    for (k=0; k<ntest; k++) {
	initVars(x, y, mat, ln, gn, sr);
	t = mvallgather(x, mat, y, gn, ln, sr, er, &t2);
	if (t < t_ag) t_ag = t;
	if (t2 < tcomm_ag) tcomm_ag = t2;
	if (verbose && gn < 80) printVec(y, ln);

	initVars(x, y, mat, ln, gn, sr);
	t = mvring(x, mat, y, gn, ln, sr, er, &t2);
	if (t < t_ring) t_ring = t;
	if (t2 < tcomm_ring) tcomm_ring = t2;
	if (verbose && gn < 80) printVec(y, ln);
    }
    if (wrank == 0) {
	double nops, r_ag, r_ring;
	printf("%d\t%d\t%.2e\t%.2e\n", gn, wsize, t_ag, t_ring);
	printf("\t%.2e\t%.2e\n", tcomm_ag, tcomm_ring);
	nops   = 2 * (double)gn * (double)gn;
	r_ag   = nops / t_ag / wsize;
	r_ring = nops / t_ring / wsize;
	printf("Perf per core\t%.2e\t%.2e\n", r_ag, r_ring);
    }

    free(x); free(y); free(mat);
    MPI_Finalize();
    return 0;
}

double mvallgather(const double *x, const double *mat, double *y,
		   int gn, int ln, int sr, int er, double *tcomm)
{
    double *gx, t, t_gather;

    gx = (double *)malloc(gn*sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    MPI_Allgather((double*)x, ln, MPI_DOUBLE, gx, ln, MPI_DOUBLE, MPI_COMM_WORLD);
    t_gather = MPI_Wtime()-t;
    mvmult(gx, mat, y, ln, gn, gn);
    t = MPI_Wtime() - t;
    *tcomm = t_gather;

    MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    free(gx);
    return t;
}

double mvring(const double *x, const double *mat, double *y,
	      int gn, int ln, int sr, int er, double *tcomm_ptr)
{
    double *x1, *x2, *xbuf, *xr, t, t2, tcomm=0.0;
    int k, wsize, wrank, nxt, prev, off;
    MPI_Request r[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    prev = wrank - 1;
    nxt  = wrank + 1;
    if (prev < 0)     prev = wsize - 1;
    if (nxt >= wsize) nxt = 0;

    x1 = (double *)malloc(ln*2*sizeof(double));
    x2 = x1 + ln;

    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();

    off  = wrank * ln;
    xbuf = (double *)x;
    xr   = x1;
    for (k=0; k<wsize-1; k++) {
	MPI_Irecv(xr, ln, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &r[0]);
	MPI_Isend(xbuf, ln, MPI_DOUBLE, nxt, 0, MPI_COMM_WORLD, &r[1]);
	mvmult(xbuf, mat+off, y, ln, ln, gn);
	t2 = MPI_Wtime();
	MPI_Waitall(2, r, MPI_STATUSES_IGNORE);
	tcomm += MPI_Wtime() - t2;
	off -= ln;  if (off < 0) off = (wsize - 1)*ln;
	xbuf = xr;
	if (xr == x1) xr = x2; else xr = x1;
    }
    /* Last time, we don't need to communicate */
    mvmult(xbuf, mat+off, y, ln, ln, gn);
    t = MPI_Wtime() - t;

    *tcomm_ptr = tcomm;
    MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    free(x1);
    return t;
}

/* Local Matrix-vector multiply.
   Matrix stored by rows
*/
void mvmult(const double  x, const double  mat,
	    double  y, int rows, int cols, int colsdim)
{
    int i, j;
    register double sum;

    for (j=0; j<rows; j++) {
	sum = y[j];
	for (i=0; i<cols; i++)
	    sum += mat[i]*x[i];
	mat += colsdim;
	y[j] = sum;
    }
}

/*
 */
void initVars(double *x, double *y, double *mat, int ln, int gn, int sr)
{
    int i, j;
    for (i=0; i<ln; i++) {
	x[i] = sr + i;
	y[i] = 0.0;
	for (j=0; j<gn; j++)
	    mat[i*gn+j] = (i+sr)*gn + j;
    }
}

void printVec(double *y, int ln)
{
    double *gy;
    int    wsize, wrank, i, gn;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    gy = (double *)malloc(ln*wsize*sizeof(double));
    MPI_Gather(y, ln, MPI_DOUBLE, gy, ln, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (wrank == 0) {
	printf("Vec:\n");
	gn = ln * wsize;
	for (i=0; i<gn; i++) {
	    printf("%e\n", gy[i]);
	}
    }
    free(gy);
}