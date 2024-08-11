// weighted algorithm
//
// Optimized main loop
// Only uses one while-loop to control the entire algorithm.
// After processing one node, the algorithm either moves to the next (linearly) if there is no
// misplaced suitor-node or it picks up the suitor node for processing.

void sweight3(int n,int *ver,int *edges,int *s,double *ws,double *weight) {


  int i,j;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;
  int count,current;

// Start of matching algorithm

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
  }
//  count = 0;

  i = 1;
  current = i;

  while (current != n+1) {
    double heaviest = ws[current];          // Find the heaviest candidate that does not have a better offer
    partner = s[current];                  // No point in trying a partner worse than the best suitor

    for(j=ver[current];j<ver[current+1];j++) { // Loop over neighbors of the current vertex 
      int y = edges[j];            // y is the neighbor of the current vertex
      if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(current,y) is the best so far, and if it is a better option for y
        heaviest = weight[j];      // Store the weight of the heaviest edge found so far
        partner = y;               // Store the name of the associated neighbor
      } // loop over neighbors
    }

    
    if (heaviest > 0) {            // True if there is a new partner
      if (s[partner] != 0) {       // True if partner already had a suitor
        next_vertex = s[partner];  // Pick up the old suitor and continue
      }
      else {
        i = i + 1;                 // Move to the next vertex
        next_vertex = i;
      }
      s[partner] = current;        // current is now the current suitor of s
      ws[partner] = heaviest;      // The weight of the edge (current,partner) 
    }
    else {
      i = i + 1;                   // Move to the next vertex
      next_vertex = i;
    }
    current = next_vertex;         // Update current vertex with the next one to work on
  } // loop over vertices

/*
   printf("Average number of iterations %f \n",(double)count/(double)n);

  int problems = 0;
  for(i=1;i<=n;i++) {
    if ((s[i] != 0) && (s[s[i]] != i))
      problems++;
  }
  if (problems > 0)
    printf("Have %d nodes with problems\n",problems);
*/
}
