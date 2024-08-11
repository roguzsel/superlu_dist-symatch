// Sequential greedy unweighted matching
// 
// Traverses the edges in their given order and use an edge if both endpoints are unused

void eseq1(int n,int m,edge *e,int *p) {

  int i,j;
  int x,y;

//  int nr_matched = 0;

// Start of matching algorithm

//  printf("Starting sequential edge based matching algorithm \n");

  for(i=1;i<=n;i++)   // Set that every node is unmatched
    p[i] = 0;

  for(i=0;i<m;i++) {                  // Loop over edges
    x = e[i].x;
    y = e[i].y;
    if ((p[x] == 0) && (p[y] == 0)) { // Check if endpoints of edge are unmatched
          p[x] = y;                   // Match x and y
          p[y] = x;                   // Match x and y
//          nr_matched += 2;
    }
  } // loop over edges

  
//  printf("Edge based greedy has matched %d vertices out of %d \n",nr_matched,n);
}
