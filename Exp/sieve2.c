#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int main( int argc, char **argv ) {

    double time_start, time_end, BW;
    struct timeval tv;
    struct timezone tz;
    int i,k=2;
    int j;
    int *prime;
    int n; //size of the Matrix.
    int approach = 0;
    /* 1. Parse arguments */
    
    int expected_argc = 1+1+1; // Expected number of arguments (executable name counts)
    if (argc != expected_argc) {
        fprintf(stderr, "Usage: ./matrix_transpose_simple N\n");
        fprintf(stderr, "ERROR: Received %d arguments, expecting %d\n", argc, expected_argc);
        exit(-1);
    }
    
    n = atoi(argv[1]);
    approach=atoi(argv[2]);
    
    fprintf(stderr,"Inputs:\n");
    fprintf(stderr,"  N=%d\n",n);
    
    // Safeguard on the allowed limit
    // Overwrite if necessary to reach higher sizes
    long memSize = sizeof(float)*n*2; // Two NxN float matrices
    long maxSize = 8L*1024*1024*1024; // 8 GiB limit (up to RAM size on the platform, to avoid swapping to hard disk)
    if (memSize > maxSize) {
        fprintf(stderr,"ERROR: trying to allocate array of size %ld bytes > maxSize = %ld bytes. Aborted.\n", memSize, maxSize);
        fprintf(stderr,"Change maxSize in the code to overwrite\n");
        exit(-1);
    }
    
    /* 2. Prepare the experiment: allocation, initialization */
    
    prime=(int*) malloc(memSize);
    
    /* initializes the prime number array. */
    
     for(i=0;i<n;i++) // 
         {
      prime[i]=i;
        }
    
    /* 3. Launch the timed experiment */

      if (approach == 0) {// MAtrix with blocks:

       gettimeofday(&tv , &tz);
        time_start =(double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
        for(i=2;i<n;i++) // Implementation of the Sieve
        {
           if(prime[i]!=0)
           { 
              for(j=2;j<n;j++)
              {
                 {
                    prime[j*prime[i]]=0;
                    if(prime[i]*j>n)
                       break;    
                 }
              }
           }
        }
        gettimeofday (&tv , &tz);
        time_end = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
        }

   
    
    // Prints the prime numbers
  
    // Raw output to stdout
    BW=sizeof(float)*n*2/(time_end-time_start)/1000000;
    printf("%d %lf %lf  \n",n,(time_end - time_start),BW);
    free(prime);
   }     