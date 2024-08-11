#include <stdio.h>
#include <stdlib.h>
#include "frame1.h"


static int compedge(const void *m1, const void *m2) {
  wed *v1 = (wed *) m1; 
  wed *v2 = (wed *) m2; 
  return (int) (v1->w < v2->w);
}

main() {

  wed *we;

  we = (wed *) malloc(2 * sizeof(wed));

  we[0].w = 0.0;
  we[0].x = 0;
  we[0].y = 1;

  we[1].w = 0.126799;
  we[1].x = 1;
  we[1].y = 2;

  printf("Before sorting \n");
  printf("0. %f \n",we[0].w);
  printf("1. %f \n",we[1].w);

  qsort(we,2,sizeof(wed),compedge);

  printf("After sorting \n");
  printf("0. %f \n",we[0].w);
  printf("1. %f \n",we[1].w);

}
