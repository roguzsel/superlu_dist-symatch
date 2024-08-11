// Round based local max algorithm.
//
// Based on the master thesis of Marcel Birn
//

void localmax(int n,int m, wed *edgel,int *p,double *ws,int *init) {

  int i,j;
  int x,y;
  int new_m;
  double w;
  wed *temp;
  wed *new_list;

  new_list = (wed *) malloc(sizeof(wed) * m);

  for(i=0;i<=n;i++) {
    p[i] = 0;        // Set that vertex i is unmatched
    ws[i] = -1;      // Set the weight of the best candidate
  }
  int r=1;
  int old_m = m;
  while (m > 0) {
  //  if (r>150) return;
  //  printf("Round number %d, remaining edges %d, removed %d \n",r,m,old_m - m);
  //  r++;
    
// First find the best candidate for each vertex

    for(i=0;i<m;i++) {
      x = edgel[i].x;
      y = edgel[i].y;
      w = edgel[i].w;

      if (w > ws[x]) { // Check if x gets a better value with this edge
        ws[x] = w;
        p[x] = y;
      }

      else if ((w == ws[x]) && (p[p[x]] != x) && (p[x]!=0)) {
        ws[x] = w;
        p[x] = y;
      }


      if (w > ws[y]) { // Check if y gets a better value with this edge
        ws[y] = w;
        p[y] = x;
      }

      else if ((w == ws[y]) && (p[p[y]] != y) && (p[y]!=0)) {
        ws[y] = w;
        p[y] = x;
      }

    }

// Go through the edges and keep those that are incident on two unmatched vertices

    new_m = 0;
    for(i=0;i<m;i++) {
      x = edgel[i].x;
      y = edgel[i].y;
      if ((p[p[x]] != x) && (p[p[y]] != y)) {  // Check if x and y are unmatched
        new_list[new_m].x = x;
        new_list[new_m].y = y;
        new_list[new_m].w = edgel[i].w;
        new_m++;
        ws[x] = -1;			// Set that x and y are unmatched
        ws[y] = -1;
        p[x] = 0;
        p[y] = 0;
      } // end if
    } // end for

    old_m = m;
    m = new_m;			// m is the number of remaining edges
    temp = edgel;		// Exchange lists
    edgel = new_list;
    new_list = temp;
  } // end while

  for(i=1;i<=n;i++) {  // Set all unmatched vertices to point to 0
    if (p[p[i]] != i)
      p[i] = 0;
  }
}
