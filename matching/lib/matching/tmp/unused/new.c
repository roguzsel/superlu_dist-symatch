#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define dim 8
#define nb 4

#define bdim  dim/nb
#define bdim2 bdim*bdim
#define bd    bdim*dim


#define n_runs 1

int main(int argc, char *argv[]) {

  double *a;    
  double *b;    
  double *c;    


  int i,j,k,i1,j1,k1;
  int ii,jf,aj,ai,xc,xr,xi,xj;
  int l;

  double mt1,mt2;
  float t_bs;
  
  printf("bdim = %d\n",bdim);
  printf("bdim2 = %d\n",bdim2);
  printf("bd = %d\n",bd);

  printf("Allocating memory \n");


  a = (double *) malloc(sizeof(double)*dim*dim);
  b = (double *) malloc(sizeof(double)*dim*dim);
  c = (double *) malloc(sizeof(double)*dim*dim);

  printf("Initializing \n");

  for(i=0;i<dim*dim;i++) {    
      a[i] = 1.0;
      b[i] = 2.0;
      c[i] = 0.0;
  }
  float flops = (float) dim* (float)dim*(2*dim-1);

  printf("Starting blocked algorithm \n");

    t_bs = -1.0;

    for(l=0;l<n_runs;l++) {
      int i1,j1,k1;
      mt1 = omp_get_wtime();


      xc = 0;
      ii = 0;
      for(i=0;i<nb;i++) {
        jf = 0;
//        printf("i = %d \n",i);
        for(j=0;j<nb;j++) {
          aj = jf;
          ai = ii;
          // printf("Block nr %d,%d \n",i,j);
          for(k=0;k<nb;k++) {
            xr = xc;
            xi = ai;
            // printf("k = %d \n",k);
            for(i1=0;i1<bdim;i1++) {
              // printf("i1 = %d \n",i1);
              xj = aj;
              for(j1=0;j1<bdim;j1++) {
                // printf("j1 = %d \n",j1);
                for(k1;k1<bdim;k1++) {               
                  // printf("k1 = %d \n",k1);
                  c[xr] += a[xi+k1] * b[xj+k1];
                } // k1
                xr++;
                xj += bdim;
              } // j1
              xi += bdim;
            } // i1
            aj += bdim2;
            ai += bdim2;
          } // k
          jf += bd;
          xc += bdim2;
        } // j
        ii += bd;
        xc = ii;
      } // i
 

      mt2 = omp_get_wtime();

      if ((t_bs < 0) || (mt2-mt1 < t_bs))
        t_bs = mt2-mt1;
    }

    printf("Matrix multiplication of two matrices of dimension %d took %f seconds\n",dim,t_bs);
    //flops = dim*dim*(2*dim-1);
    printf("Flop count is %d \n",flops);
    printf("MFlops rate is %lf \n",flops/t_bs/1000000);

    printf("%lf \n",c[0]);


}
