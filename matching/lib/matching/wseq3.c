// Sequential greedy matching
//
// Pick the heaviest available incident edge for each vertex 
// Vertices are processed as given by order[] (ie. random)
//
// This code has been optimized with prefetching directives
//
// Note that the prefetching can cause problems, most likely because at the end
// of processing where it is stretching into vertices that don't exist.

void wseq3(int n,int *ver,int *edges,int *p,double *weight, int *order ) {

  int i,j,k;
  int nr_union = 0;
  int nr_matched = 0;
  double heaviest;
  int partner;
  int y;

// Start of matching algorithm


  for(i=1;i<=n;i++)   // Set that every node is unmatched
    p[i] = 0;

  for(k=1;k<n;k+=2) {                 // Loop over vertices, no need to process last vertex

    __builtin_prefetch (&order[k+4], 0, 3);
    __builtin_prefetch (&p[order[k+2]], 1, 0);
    __builtin_prefetch (&ver[order[k+2]], 0, 3);
    __builtin_prefetch (&edges[ver[order[k+2]]], 0, 3);
    __builtin_prefetch (&weight[ver[order[k+2]]], 0, 3);

    __builtin_prefetch (&p[order[k+3]], 1, 0);
    __builtin_prefetch (&ver[order[k+3]], 0, 3);
    __builtin_prefetch (&edges[ver[order[k+3]]], 0, 3);
    __builtin_prefetch (&weight[ver[order[k+3]]], 0, 3);

    i = order[k];
      
    if (p[i] == 0) {                 // Check if this vertex has been matched
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
        // nr_matched += 2;
      }
    }

    i = order[k+1];
    if (p[i] == 0) {                 // Check if this vertex has been matched
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
        // nr_matched += 2;
      }
    }

  } // loop over vertices

  // printf("Have matched %d vertices out of %d \n",nr_matched,n);
}
