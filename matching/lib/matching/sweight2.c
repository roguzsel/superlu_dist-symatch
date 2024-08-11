// weighted algorithm
// This code is not in use

void sweight(int n,int *ver,int *edges,int *s,float *ws,float *weight) {


  int i,j;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;
  int count;

// Start of matching algorithm

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
  }
//  count = 0;

  for(x=1;x<=n;x++) {                // Loop over vertices
    int done = FALSE;
    i = x;                           // Looking for a matching partner for node i
    while (!done) {
//      count++;

      float heaviest = 0.0;          // Find the heaviest candidate that does not have a better offer

      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int y = edges[j];            // y is the current neighbor of i
        if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y
          heaviest = weight[j];      // Store the weight of the heaviest edge found so far
          partner = y;               // Store the name of the associated neighbor
        }
      } // loop over neighbors


      if (heaviest > 0) {            // True if there is a partner
        if (s[partner] != 0) {       // True if partner already had a suitor
          next_vertex = s[partner];  // Pick up the old suitor and continue
          s[partner] = i;              // i is now the current suitor of s
          i = next_vertex;
          done = FALSE;
        }
        else {
          s[partner] = i;              // i is now the current suitor of s
          done = TRUE;                 // No neighbor to match with
        }
        ws[partner] = heaviest;      // The weight of the edge (i,partner) 
      }
      else
        done = TRUE;                 // No neighbor to match with

    } // while not done
  } // loop over vertices

//  printf("Average number of iterations %f \n",(float)count/(float)n);
}
