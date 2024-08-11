// Parallel Suitor algorithm using a queue
//
// "n" is the number of vertices
// The graph is stored in "ver" and "edges" using compressed neighbor lists. "ver[i]" points into where in "edges" the neighbors of vertex i
// are stored. Note that each edge is stored twice (in both directions). Numbering of the vertices starts from 1!
// "s" is an array of length n+1 to store the suitor value of each vertex (note that position 0 is not used).
// "q1 and q2" are arrays of length n+1, used for storing the vertices that must be worked on
// "ws" stores the weight of each suitor, indexed just like "s".
// "n_locks" are OpenMP locks, one for each vertex. Used when setting the suitor value, and also when testing if a vertex can beat an
// old suitor with the same weight (ties are broken based on vertex ids).
// "weight" stores the weights of the edges, ordered just like "edges".

void pqueue1(int n,int *ver,int *edges,int *s,int *q1,int *q2,double *ws,omp_lock_t *nlocks,double *weight) {


  int i,j,skip;
  int partner;     // Prospective candidate to match with   
  int next_vertex;

  int my_id = omp_get_thread_num();
  int threads = omp_get_num_threads();

  int done,x;

// Start of matching algorithm
// First set the suitor values and corresponding weights to all 0 (i.e. unused).
// Also initialize the locks.

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {   
    s[i] = 0;                  // Set that no node is trying to match with i 
    ws[i] = 0.0;               // Set the current weight of best suitor edge to 0
    omp_init_lock(&nlocks[i]); // Initialize locks
    q1[i] = i;                 // put each node in the queue
  }

// Set up the intervals for each thread. Only static load balance.
  int start,stop;

  int chunk = n/threads;
  int rest = n%threads;
  if (my_id < rest) {
    start = 1+(chunk+1)*my_id;
    stop = start + chunk + 1;
  }
  else {
    start = 1+rest + my_id*chunk;
    stop = start + chunk;
  }

#pragma omp barrier

// Now start the real matching algorithm.
// The following loop is run once for each vertex. It is the main parallel loop.
//
  
  int count = 0;

  while (start < stop) {              // The current queue stretches from start to stop - 1
    int newfirst = start;             // Starting point of the new queue
    int first = start;                // Starting point of the current queue
    while (first < stop) {
      i = q1[first++];                // Get first element, move queue pointer one position forward 

      double heaviest = -1.0;  // All weights are positive. "heaviest" contains the weight of the best one so far.

      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int y = edges[j];            // y is the current neighbor of i
      //printf("%d is processing %d and considering %d \n",my_id,i,y);

        if ((weight[j] < heaviest) || (weight[j] < ws[y])) // If the weight of the current edge is lighter than the current best or
          continue;                                        // if we cannot offer y something better than its current suitor, move on.

        if ((weight[j] == heaviest) && (y < partner))      // When equal weight, give priority to highest index
          continue;

        if (weight[j] == ws[y])  {                         // Testing if we beat old suitor with the same weight.
          skip = false;
          omp_set_lock(&nlocks[y]);    // Locking y, needed to ensure that s[y] and ws[y] were set using the same vertex. 
            if ((weight[j] == ws[y]) && (i < s[y]))        // Must have higher index to beat previous suitor when equal weight
              skip = true;
          omp_unset_lock(&nlocks[y]); // Unlocking y
          if (skip)                                        // If index was lower then ignore this edge
            continue;
        }

        heaviest = weight[j];      // Store the weight of the best so far
        partner = y;               // Also store the vertex id.
       
      } // loop over neighbors

// Done finding a candidate. If heaviest > 0 then we try to set i as a suitor of partner.

      if (heaviest > 0) {
        omp_set_lock(&nlocks[partner]);    // Locking partner to ensure that only one vertex tries to access s[] and ws[].
// Test that we still beat the old suitor, this might have changed since last time we checked.
        if ((heaviest > ws[partner]) || ((heaviest == ws[partner]) && (i>s[partner]))) {
      
          if (s[partner] != 0) {       // True if partner already had a suitor
            q2[newfirst++] = s[partner]; // push old partner on stack, increase stack size 
            count++;
          }

          s[partner] = i;              // i is now the suitor of partner
          ws[partner] = heaviest;      // The weight of the edge (current,partner) 

          omp_unset_lock(&nlocks[partner]); // Unlocking partner
        }
        else {   // can no longer use the result for node i, must find new partner for it, continue while loop.
          omp_unset_lock(&nlocks[partner]); // Unlocking partner
          q2[newfirst++] = i;               // put i back in the queue, increase queue length
        }
      } // end if heaviest > 0
    } // end loop over queue-length
    stop = newfirst;                   // Switch the two queues and memorize stopping position of new queue
    int *t = q1;
    q1 = q2;
    q2 = t;
  }  // while the new queue is not empty
#pragma omp barrier
}
