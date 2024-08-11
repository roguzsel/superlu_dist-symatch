// Sequential greedy matching
// Goes through the vertices and picks the first unmatched neighbor to match with

void seq1(int n,int *ver,int *edges,int *p) {

  int i,j;
  int nr_union = 0;
// int nr_matched = 0;

// Start of matching algorithm

//  printf("Starting sequential matching algorithm \n");

  for(i=1;i<=n;i++)   // Set that every node is unmatched
    p[i] = 0;

  for(i=1;i<=n;i++) {                // Loop over vertices
    if (p[i] == 0) {                 // Check if this vertex has been matched
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        if (p[edges[j]] == 0) {      // Check if vertex is unmatched
          p[i] = edges[j];           // Match i and e[j]
          p[edges[j]] = i;
//          nr_matched += 2;
          break;                     // Continue to next vertex
        }
      }
    }
  } // loop over vertices

  
//  printf("Have matched %d vertices out of %d \n",nr_matched,n);
}
