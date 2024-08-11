// weighted algorithm
//

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
  count = 0;

  for(x=1;x<=n;x++) {                // Loop over vertices
    int done = FALSE;
    i = x;                           // Looking for a matching partner for node i
    while (!done) {
      count++;

      float heaviest = ws[i];          // Find the heaviest candidate that does not have a better offer
      partner = s[i];                  // No point in trying a partner worse than the best suitor
      
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int y = edges[j];            // y is the current neighbor of i

//        if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y

        if ((weight[j] < heaviest) || (weight[j] < ws[y]))
          continue;

        if ((weight[j] == heaviest) && (y < partner))
          continue;

        if ((weight[j] == ws[y]) && (i < s[y]))
          continue;

          heaviest = weight[j];      // Store the weight of the heaviest edge found so far
          partner = y;               // Store the name of the associated neighbor
//        }
      } // loop over neighbors

      done = TRUE;                   // Check if we found someone to match with

      if (heaviest > 0) {            // True if there is a partner
        if (s[partner] != 0) {       // True if partner already had a suitor
          next_vertex = s[partner];  // Pick up the old suitor and continue
          done = FALSE;
        }
        s[partner] = i;              // i is now the current suitor of s
        ws[partner] = heaviest;      // The weight of the edge (i,partner) 
      }
      if (!done) {
        i = next_vertex;
      }
    } // while not done
  } // loop over vertices

  printf("Average number of iterations %f \n",(float)count/(float)n);
}
