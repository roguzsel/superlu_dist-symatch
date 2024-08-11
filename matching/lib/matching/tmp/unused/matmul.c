/*        Basic matrix multiplication algorithm. Computes C = A*B
          where A,B, and C are dim x dim matrices. The program prints
          out the time taken, the number of flops executed, and the
          number of flops per second.

          Compile with gcc mult.m -fopenmp -O3
          Run with ./a.out
*/
 
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define dim 1024    // dimension of matrix
#define n_runs 1    // number of times to run the main algorithm

int main(int argc, char *argv[]) {

  double **a;    
  double **b;    
  double **c;    

  double *aptr;
  double *bptr;
  double *cptr;

  int i,j,k;
  int l;

  double mt1,mt2;
  float t_bs;

  printf("Allocating memory \n");

  a =  malloc(sizeof(double)*dim);
  b =  malloc(sizeof(double)*dim);
  c =  malloc(sizeof(double)*dim);

  aptr = (double *) malloc(sizeof(double)*dim*dim);
  bptr = (double *) malloc(sizeof(double)*dim*dim);
  cptr = (double *) malloc(sizeof(double)*dim*dim);

  for(i=0;i<dim;i++) {   // Setting pointers
    a[i] = aptr + i*dim; 
    b[i] = bptr + i*dim;
    c[i] = cptr + i*dim;
  }

  printf("Initializing \n");

  for(i=0;i<dim;i++) {    // Initializing data
    for(j=0;j<dim;j++) {     
      a[i][j] = 1.0;
      b[i][j] = 2.0;
      c[i][j] = 0.0;
    }
  }

  float flops = (float) dim* (float)dim*(2*dim-1);  // Calculate the number of flops

  printf("Starting algorithm \n");

    t_bs = -1.0;   // Set best running time to negative value to assure that we pick up the first run

    for(l=0;l<n_runs;l++) {   // Run the algorithm n_runs times

      mt1 = omp_get_wtime();  // Start of timing
 
      for(i=0;i<dim;i++) {    
        for(j=0;j<dim;j++) {     
          for(k=0;k<dim;k++) {     
            c[i][j] += a[i][k] * b[k][j];
          }
        }
      }

      mt2 = omp_get_wtime(); // End of timing

      if ((t_bs < 0) || (mt2-mt1 < t_bs))  // Storing the best run, could also have taken average
        t_bs = mt2-mt1;
    }

    printf("Matrix multiplication of two matrices of dimension %d took %f seconds\n",dim,t_bs);
    printf("Flop count is %.0f \n",flops);
    printf("MFlop rate is %lf \n",flops/t_bs/1000000);

    printf("%lf \n",c[0][0]); // Print out something to fool the compiler

}
