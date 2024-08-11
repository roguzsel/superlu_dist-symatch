// Sequential greedy matching
//
// Pick the heaviest available incident edge for each vertex 
// Vertices are processed in order

void wseq1(int n,int *ver,int *edges,int *p,double *weight ) {

  int i,j;
  int nr_union = 0;
  int nr_matched = 0;
  double heaviest;
  int partner;
  int y;

// Start of matching algorithm

//  printf("Starting sequential matching algorithm \n");

  for(i=1;i<=n;i++)   // Set that every node is unmatched
    p[i] = 0;

  for(i=1;i<=n;i++) {                 // Loop over vertices
    if (p[i] == 0) {                  // Check if this vertex has been matched
      heaviest = 0.0;
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        y = edges[j];

        if ((p[y] == 0) && (weight[j]>heaviest)) {     // Check if vertex is unmatched and if it gives a new heaviest edge
          heaviest = weight[j];
          partner = y;
        }
      }

      if (heaviest > 0.0) {
        p[i] = partner;           // Match i and partner
        p[partner] = i;
  //      nr_matched += 2;
      }

    }
  } // loop over vertices
  
//  printf("Have matched %d vertices out of %d \n",nr_matched,n);
}
