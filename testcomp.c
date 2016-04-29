#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


void rank_printf(int rank, const char* fmt, ...)
{
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if (myRank == rank) {
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
    }
}

void pingpong_blocking(int i, int j, int m, int rank, char* buffer)
{
    // printf("In func, my rank is %d, i=%d, j=%d\n", rank, i, j);
    if (rank == i) {
        MPI_Send(buffer, m, MPI_CHAR, j, 0, MPI_COMM_WORLD);
        MPI_Recv(buffer, m, MPI_CHAR, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // printf("Sent from process %d\n", i);
    }
    if (rank == j) {
        MPI_Recv(buffer, m, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer, m, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        // printf("Received in process %d\n", j);
    }
}

void pingpong_nonblocking(int i, int j, int m, int rank, char* buffer)
{
    // printf("In func, my rank is %d, i=%d, j=%d\n", rank, i, j);
    if (rank == i) {
        MPI_Request req;
        MPI_Isend(buffer, m, MPI_CHAR, j, 0, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);

        MPI_Irecv(buffer, m, MPI_CHAR, j, 0, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        // printf("Sent from process %d\n", i);
    }
    if (rank == j) {
        MPI_Request req;
        MPI_Irecv(buffer, m, MPI_CHAR, i, 0, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);

        MPI_Isend(buffer, m, MPI_CHAR, i, 0, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        // printf("Received in process %d\n", j);
    }
}

void head_to_head(int i, int j, int m, int rank, char* msg_out, char* msg_in)
{
    if ((rank != i) && (rank != j)) return;

    int other = (rank == i ? j : i);

    MPI_Request req[2];
    MPI_Isend(msg_out, m, MPI_CHAR, other, 0, MPI_COMM_WORLD, &req[0]);
    MPI_Irecv(msg_in, m, MPI_CHAR, other, 0, MPI_COMM_WORLD, &req[1]);
    MPI_Waitall(2, req, MPI_STATUS_IGNORE);
}

void print_header(void)
{
    rank_printf(0, "test\ti\tj\tm\ttime\n");
}

void run_tests(int i, int j, int rank, int *mvals, int mcount, int niters, int nPerIter)
{
    double start, end;

    for (int mIter = 0; mIter < mcount; mIter++) {
        int m = mvals[mIter];

        char* buffer1 = malloc(m * sizeof(char));
        char* buffer2 = malloc(m * sizeof(char));

        if (rank == i || rank == j) {
            double minTime = 1e9;
            for (int iter = 0; iter < niters; iter++) {
                start = MPI_Wtime();
                for (int run = 0; run < nPerIter; run++) {
                    pingpong_blocking(i, j, m, rank, buffer1);
                }
                end = MPI_Wtime();
                double timeTaken = (end - start) / nPerIter;
                if (timeTaken < minTime) minTime = timeTaken;
            }
            rank_printf(i, "ppnb\t%d\t%d\t%d\t%0.4le\n", i, j, m, minTime);
        }

        MPI_Barrier(MPI_COMM_WORLD);  // ----------------------------------------------------------------------

        if (rank == i || rank == j) {
            double minTime = 1e9;
            for (int iter = 0; iter < niters; iter++) {
                start = MPI_Wtime();
                for (int run = 0; run < nPerIter; run++) {
                    pingpong_nonblocking(i, j, m, rank, buffer1);
                }
                end = MPI_Wtime();
                double timeTaken = (end - start) / nPerIter;
                if (timeTaken < minTime) minTime = timeTaken;
            }
            rank_printf(i, "ppb\t%d\t%d\t%d\t%0.4le\n", i, j, m, minTime);
        }

        MPI_Barrier(MPI_COMM_WORLD);  // ----------------------------------------------------------------------

        if (rank == i || rank == j) {
            double minTime = 1e9;
            for (int iter = 0; iter < niters; iter++) {
                start = MPI_Wtime();
                for (int run = 0; run < nPerIter; run++) {
                    head_to_head(i, j, m, rank, buffer1, buffer2);
                }
                end = MPI_Wtime();
                double timeTaken = (end - start) / nPerIter;
                if (timeTaken < minTime) minTime = timeTaken;
            }
            rank_printf(i, "h2h\t%d\t%d\t%d\t%0.4le\n", i, j, m, minTime);
        }

        MPI_Barrier(MPI_COMM_WORLD);  // ----------------------------------------------------------------------

        free(buffer1);
        free(buffer2);
    }
}

int main(int argc, char** argv)
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rank_printf(0, "Running with %d MPI processes\n", size);

    if (rank == 0) {
        double tickSize = MPI_Wtick();
        printf("MPI wall clock resolution: %0.4le seconds\n", tickSize);
    }

    int mvals[18];
    mvals[0] = 0;
    for (int k = 2; k <= 18; k++) {
        mvals[k-1] = (1 << k);
    }

    rank_printf(0, "\n");

    const int numIters = 5;
    const int numPerIter = 1000;

    print_header();

    run_tests(0, 1, rank, mvals, 18, numIters, numPerIter);
    run_tests(0, size-1, rank, mvals, 18, numIters, numPerIter);
    run_tests(2, 3, rank, mvals, 18, numIters, numPerIter);

    MPI_Finalize();
    return 0;
}