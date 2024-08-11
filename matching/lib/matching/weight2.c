// Parallel Suitor algorithm
//
// "n" is the number of vertices
// The graph is stored in "ver" and "edges" using compressed neighbor lists. "ver[i]" points into where in "edges" the neighbors of vertex i
// are stored. Note that each edge is stored twice (in both directions). Numbering of the vertices starts from 1!
// "s" is an array of length n+1 to store the suitor value of each vertex (note that position 0 is not used).
// "ws" stores the weight of each suitor, indexed just like "s".
// "n_locks" are OpenMP locks, one for each vertex. Used when setting the suitor value, and also when testing if a vertex can beat an
// old suitor with the same weight (ties are broken based on vertex ids).
// "weight" stores the weights of the edges, ordered just like "edges".

void weighted2(int n,int *ver,int *edges,int *s,double *ws,omp_lock_t *nlocks,double *weight) {


  int i,j,skip;
  int partner;     // Prospective candidate to match with   
  int next_vertex;

//  int my_id = omp_get_thread_num();
//  int threads = omp_get_num_threads();

  int done,x;

// Start of matching algorithm
// First set the suitor values and corresponding weights to all 0 (i.e. unused).
// Also initialize the locks.

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {   
    s[i] = 0;                  // Set that no node is trying to match with i 
    ws[i] = 0.0;               // Set the current weight of best suitor edge to 0
    omp_init_lock(&nlocks[i]); // Initialize locks
  }

#pragma omp barrier

// Now start the real matching algorithm.
// The following loop is run once for each vertex. It is the main parallel loop.

// #pragma omp for schedule(dynamic,700) private(i,j,done,x,partner,next_vertex)
#pragma omp for schedule(static) private(i,j,done,x,partner,next_vertex)
  for(x=1;x<=n;x++) {                // Loop over vertices
    i = x;                           // i is the vertex that we are trying to make the suitor of one of its neighbors
    done = false;                    // done is used to ensure that we have finished processing any old suitor that might get dislodged.

    while (!done) {

/*  // This was used for testing if starting the search from the current suitor was faster. 
      omp_set_lock(&nlocks[i]);    // Locking myself
      double heaviest = ws[i];   // Start search from current suitor
      partner = s[i];
      omp_unset_lock(&nlocks[i]); // Unlocking myself
*/

      double heaviest = -1.0;  // All weights are positive. "heaviest" contains the weight of the best one so far.

      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int y = edges[j];            // y is the current neighbor of i
/*
        if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y
*/


        if ((weight[j] < heaviest) || (weight[j] < ws[y])) // If the weight of the current edge is lighter than the current best or
          continue;                                        // if we cannot offer y something better than its current suitor, move on.

        if ((weight[j] == heaviest) && (y < partner))      // When equal weight, give priority to highest index
          continue;

/*
        if ((weight[j] == ws[y]) && (i < s[y]))            // Must have higher index to beat previous suitor when equal weight
          continue;
*/
        if (weight[j] == ws[y])  {                         // Testing if we beat old suitor with the same weight.
          skip = false;
          omp_set_lock(&nlocks[y]);    // Locking y, needed to ensure that s[y] and ws[y] were set using the same vertex. 
            if ((weight[j] == ws[y]) && (i < s[y]))        // Must have higher index to beat previous suitor when equal weight
              skip = true;
          omp_unset_lock(&nlocks[y]); // Unlocking y
          if (skip)                                        // If index was lower then ignore this edge
            continue;
        }

// Check if w(i,y) is the best so far, and if it is a better option for y
/*
        if (((weight[j] > heaviest) || ((weight[j] == heaviest) && (y > partner))) &&
            ((weight[j] > ws[y]) || ((weight[j] == ws[y]) && (i > s[y])))) {
*/
          heaviest = weight[j];      // Store the weight of the best so far
          partner = y;               // Also store the vertex id.
//       }
       
      } // loop over neighbors

// Done finding a candidate. If heaviest > 0 then we try to set i as a suitor of partner.

      done = true;
      if (heaviest > 0) {
        omp_set_lock(&nlocks[partner]);    // Locking partner to ensure that only one vertex tries to access s[] and ws[].

// Test that we still beat the old suitor, this might have changed since last time we checked.
        if ((heaviest > ws[partner]) || ((heaviest == ws[partner]) && (i>s[partner]))) {
        // if (heaviest >= ws[partner]) { // Must have >= because i might start from its best suitor
          next_vertex = s[partner]; // Store the old suitor
          s[partner] = i;           // i is now the current suitor of s
          ws[partner] = heaviest;   // The weight of the edge (i,partner) 
          omp_unset_lock(&nlocks[partner]); // Unlocking partner

          if (next_vertex != 0) {   // If there was an old suitor then we must process it
            i = next_vertex;        // i now holds the id of the old suitor
            done = false;           // This ensures that the while loop will continue 
          }
        }
        else {   // can no longer use the result for node i, must find new partner for it, continue while loop.
          omp_unset_lock(&nlocks[partner]); // Unlocking partner
          done = false;
        }
      }
    }  // while not done
  } // loop over vertices

}
