// weighted algorithm
//
// Suitor based algorithm, asuming sorted edge lists
// Parallel version

void psweight1(int n,int *ver,int *edges,int *s,double *ws,double *weight,int *next,omp_lock_t *nlocks) {


  int i,j,e,x,y;
  int done, next_vertex;


// Start of matching algorithm

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
    next[i] = ver[i];   // Set where in edge list to start searching for a partner 
    omp_init_lock(&nlocks[i]); // Initialize locks
  }
#pragma omp barrier

// #pragma omp for schedule(dynamic,700) private(x)
#pragma omp for schedule(static) private(x)
  for(x=1;x<=n;x++) {       // Loop over all vertices
    done = false;
    i = x;
    while (!done) {         // Trying to find a (possibly new) partner for i
      y = next[i];          // y points to the next possible neighbor (in the edge list of i) that i can match with
      e = edges[y];         // Get the id of the neighbor that we will try to match with


// Stop if there are no more edges or if we found a candidate
/*      while ((y < ver[i+1]) && ((weight[y] < ws[e]) || ((weight[y] == ws[e]) && (i > s[e])))) {  
        y++;              // Move to next position in edge list
        e = edges[y];     // Get id of neighbor
      }
*/
      while ((y < ver[i+1]) && (weight[y] < ws[e])) {
        y++;              // Move to next position in edge list
        e = edges[y];     // Get id of neighbor
      }

      done = true;
      if (y < ver[i+1]) {            // Test if we found a possible partner
        next[i] = y+1;               // Set where to search from next time i needs a partner
        omp_set_lock(&nlocks[e]);    // Locking e 
        if ((weight[y] > ws[e]) || ((weight[y] == ws[e]) && (i > s[e]))) {
          if (s[e] != 0) {             // True if e already had a suitor
            next_vertex = s[e];        // Pick up the old suitor
            done = false;
          }
          s[e] = i;                  // i is now the current suitor of e
          ws[e] = weight[y];           // Store the weight of (i,e) as the weight given by the current suitor for e
        }
        else {   // can no longer use the result for node i, must find new partner
          done = false;
          next_vertex = i;
        }
        omp_unset_lock(&nlocks[e]);  // Unlocking e 
      } // if y
      if (!done) {  // Continue with the next vertex
        i = next_vertex;
      }
    } // while not done
  } // loop over vertices

}
