// weighted algorithm
//
// Assuming sorted edge lists!

// Optimized main loop
// Only uses one while-loop to control the entire algorithm.
// After processing one node, the algorithm either moves to the next (linearly) if there is no
// misplaced suitor-node or it picks up the suitor node for processing.

void sweight7(int n,int *ver,int *edges,int *s,double *ws,double *weight,int *next) {

  int i,y,e;
  int next_vertex;
  int current;

// Start of matching algorithm

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
    next[i] = ver[i];   // Set where to start searching for next candidate
  }
//  count = 0;

  i = 1;
  current = i;

  while (current != n+1) {

    y = next[current];    // y points to the next possible neighbor (in the edge list of current) that current can match with
    e = edges[y];         // Get the id of the neighbor that we will try to match with

// Stop if there are no more edges or if we found a candidate
    while ((y < ver[current+1]) && ((weight[y] < ws[e]) || ((weight[y] == ws[e]) && (current > s[e])))) {
      y++;              // Move to next position in edge list
      e = edges[y];     // Get id of neighbor
    }

    if (y < ver[current+1]) {      // Test if we found a possible partner

      next[current] = y+1;         // Set where to search from next time current needs a partner

      if (s[e] != 0) {             // True if e already had a suitor
        next_vertex = s[e];        // Pick up the old suitor
      }
      else {
        i = i + 1;                 // Move to the next vertex
        next_vertex = i;
      }
      s[e] = current;              // current is now the suitor of e
      ws[e] = weight[y];           // Store the weight of (i,e) as the weight given by the current suitor for e
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
