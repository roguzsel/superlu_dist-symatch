// Parallel Suitor algorithm
//
// This version uses only value to keep track of the suitor of a vertex. The suitor value
// is then a pointer into the edge array where both the id and the weight of the suitor can be found.

void weighted2_os(int n,int *ver,int *edges,int *s,double *ws,omp_lock_t *nlocks,double *weight,int *where) {


  int i,j,skip;
  int partner;     // Prospective candidate to match with   
  int next_vertex;
  int place;

//  int my_id = omp_get_thread_num();
//  int threads = omp_get_num_threads();

  int done,x;

// Start of matching algorithm

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {   
    s[i] = ver[n+1];           // Set that no node is trying to match with i 
    // ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
    omp_init_lock(&nlocks[i]); // Initialize locks
  }

  weight[ver[n+1]] = 0.0;
  edges[ver[n+1]] = 0;

#pragma omp barrier

// #pragma omp for schedule(dynamic,700) private(i,j,done,x,partner,next_vertex)
#pragma omp for schedule(static) private(i,j,done,x,partner,next_vertex)
  for(x=1;x<=n;x++) {                // Loop over vertices
    i = x; 
    done = false;
    while (!done) {

      int who = s[i];
      double heaviest = -1.0;

      if (who != ver[n+1]) {
        heaviest = weight[who];   // Start search from current suitor
        partner = edges[who];
        place = where[who];  // where is node i placed in the neighbor list of partner
      }

      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int y = edges[j];            // y is the current neighbor of i
        who = s[y];
        if ((weight[j] < heaviest) || (weight[j] < weight[s[y]])) // If the weight of the current edge is lighter than the current best or
          continue;                                        // if we cannot offer y something better, move on.

        if ((weight[j] == heaviest) && (y <= partner))      // When equal weight, give priority to highest index
          continue;

        if ((weight[j] == weight[who]) && (i < edges[who]))            // Must have higher index to beat previous suitor when equal weight
          continue;



// Check if w(i,y) is the best so far, and if it is a better option for y
/*
        if (((weight[j] > heaviest) || ((weight[j] == heaviest) && (y > partner))) &&
            ((weight[j] > ws[y]) || ((weight[j] == ws[y]) && (i > s[y])))) {
*/
          heaviest = weight[j];      // Store the best so far
          partner = y;
          place = where[j];          // position of i in the neighbor list of y
//          printf("%d med vekt %lf er beste partner for %d \n",partner,heaviest,i);
//       }
       
      } // loop over neighbors

      done = true;
      if (heaviest > 0) {
        omp_set_lock(&nlocks[partner]);    // Locking partner

        // if ((heaviest > ws[partner]) || ((heaviest == ws[partner]) && (i>s[partner]))) {
        if (heaviest >= weight[s[partner]]) { // Must have >= because i might start from its best suitor
          if (edges[s[partner]] != 0) {
            next_vertex = edges[s[partner]];
            done = false;
          }
          s[partner] = place;
        }
        else {   // can no longer use the result for node i, must find new partner
          done = false;
          next_vertex = i;
        }
        omp_unset_lock(&nlocks[partner]); // Unlocking partner
      }
      if (!done) {  // Continue with the next vertex
        i = next_vertex;
      }
    }  // while not done
  } // loop over vertices

}
