#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define dim 1024
#define BlockSize 2
#define min(a,b) (((a) < (b)) ? (a) : (b))

#define n_runs 1

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

  for(i=0;i<dim;i++) {
    a[i] = aptr + i*dim; 
    b[i] = bptr + i*dim;
    c[i] = cptr + i*dim;
  }

  printf("Initializing \n");

  for(i=0;i<dim;i++) {    
    for(j=0;j<dim;j++) {     
      a[i][j] = 1.0;
      b[i][j] = 2.0;
      c[i][j] = 0.0;
    }
  }
  float flops = (float) dim* (float)dim*(2*dim-1);
/*
  printf("Starting algorithm \n");

    t_bs = -1.0;

    for(l=0;l<n_runs;l++) {

      mt1 = omp_get_wtime(); 
 
      for(i=0;i<dim;i++) {    
        for(j=0;j<dim;j++) {     
          for(k=0;k<dim;k++) {     
            c[i][j] += a[i][k] * b[k][j];
          }
        }
      }

      mt2 = omp_get_wtime(); 

      if ((t_bs < 0) || (mt2-mt1 < t_bs))
        t_bs = mt2-mt1;
    }

    printf("Matrix multiplication of two matrices of dimension %d took %f seconds\n",dim,t_bs);
    printf("Flop count is %d \n",flops);
    printf("MFlops rate is %lf \n",flops/t_bs/1000000);

    printf("%lf \n",c[0][0]);

*/
  printf("Starting blocked algorithm \n");

    t_bs = -1.0;

    for(l=0;l<n_runs;l++) {
      int i1,j1,k1;
      mt1 = omp_get_wtime();


      c[0][0]=0;
      for( i1=0;i1<(dim/BlockSize);++i1) {
        for(j1=0;j1<(dim/BlockSize);++j1) {
          for(k1=0;k1<(dim/BlockSize);++k1) {
            for(i=i1=0;i<min(i1+BlockSize-1);++i) {
              for(j=j1=0;j<min(j1+BlockSize-1);++j) {
                for(k=k1;k<min(k1+BlockSize-1);++k) {               
                  c[i][j] += a[i][k] * b[k][j];
                }
              }
            }
          }
        }
      }
 

      mt2 = omp_get_wtime();

      if ((t_bs < 0) || (mt2-mt1 < t_bs))
        t_bs = mt2-mt1;
    }

    printf("Matrix multiplication of two matrices of dimension %d took %f seconds\n",dim,t_bs);
    //flops = dim*dim*(2*dim-1);
    printf("Flop count is %d \n",flops);
    printf("MFlops rate is %lf \n",flops/t_bs/1000000);

    printf("%lf \n",c[0][0]);


}
