// weighted algorithm
//
// McVittie-Wilson algorithm using a stack

void mcw1(int n,int *ver,int *edges,int *s,double *ws,int *stack,double *weight) {


  int i,j;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;
  int count,current;

// Start of matching algorithm

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
    stack[i] = i;       // push each node on the stack
  }
  // count = 0;
 
  int top = 1;

  while (top != n+1) {
    current = stack[top++];                // Get top element, move stack pointer
    double heaviest = ws[current];         // Find the heaviest candidate that does not have a better offer
    partner = s[current];                  // No point in trying a partner worse than the best suitor

    for(j=ver[current];j<ver[current+1];j++) { // Loop over neighbors of the current vertex 
      int y = edges[j];            // y is the neighbor of the current vertex
// Check if w(current,y) is the best so far, and if it is a better option for y
      if ((weight[j]>heaviest) && (ws[y]<weight[j])) {
        heaviest = weight[j];      // Store the weight of the heaviest edge found so far
        partner = y;               // Store the name of the associated neighbor
      } // loop over neighbors
    }
    
    if (heaviest > 0) {            // True if there is a new partner
      if (s[partner] != 0) {       // True if partner already had a suitor
        stack[--top] = s[partner]; // push old partner on stack, increase stack size 
        // count++;
      }
      
      s[partner] = current;        // current is now the suitor of partner
      ws[partner] = heaviest;      // The weight of the edge (current,partner) 
    }
    
  } // loop over vertices

  // printf("Number of insertions is %d\n",count);
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
