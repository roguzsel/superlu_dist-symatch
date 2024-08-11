// weighted algorithm
//
// Suitor based algorithm, asuming sorted edge lists

void sweight1(int n,int *ver,int *edges,int *s,double *ws,double *weight,int *next) {


  int i,j,e,x,y;
  int done, next_vertex;

// Start of matching algorithm

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
    next[i] = ver[i];   // Set where in edge list to start searching for a partner 
  }

  int count = 0;

  for(x=1;x<=n;x++) {       // Loop over all vertices
    done = false;
    i = x;


    while (!done) {         // Trying to find a (possibly new) partner for i

      count++;
      y = next[i];          // y points to the next possible neighbor (in the edge list of i) that i can match with
      e = edges[y];         // Get the id of the neighbor that we will try to match with

// Stop if there are no more edges or if we found a candidate
      while ((y < ver[i+1]) && ((weight[y] < ws[e]) || ((weight[y] == ws[e]) && (i < s[e])))) {  
        y++;              // Move to the next position in the edge list
        e = edges[y];     // Get the id of neighbor
      }

      done = true;
      if (y < ver[i+1]) {            // Test if we found a possible partner
        next[i] = y+1;               // Set where to search from next time i needs a partner
        ws[e] = weight[y];           // Store the weight of (i,e) as the weight given by the current suitor for e

        if (s[e] != 0) {             // True if e already had a suitor
          next_vertex = s[e];        // Pick up the old suitor
          s[e] = i;                  // i is now the current suitor of e
          i = next_vertex;           // Continue and try to match the previous suitor of e
          done = false;
        }
        else {
          s[e] = i;                  // Set i as the current suitor of e
        }
      }
    } // while not done
  } // loop over vertices
  // printf("Number of vertices considered %d \n",count);
}
