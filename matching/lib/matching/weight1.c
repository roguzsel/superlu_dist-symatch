// weighted algorithm
//
//  This implementation uses a local stack on each processor to store vertices that must be processed. The 
//  processing happens in two stages. First, the vertices in the global list are processed. Any unmatched
//  vertex is then put into the stack which is handled in the second loop.


void weighted(int n,int *ver,int *edges,int *s,double *ws,int *left,omp_lock_t *nlocks,double *weight) {


  int i,j;
  int partner;     // Prospective candidate to match with   
  int nr_union = 0;
  int nr_matched = 0;

  int my_id = omp_get_thread_num();
  int threads = omp_get_num_threads();

  int start = my_id*n/threads;   // Calculate where to start storing elements in "left"
  int num_left = 0;

// Start of matching algorithm

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
    omp_init_lock(&nlocks[i]); // Initialize locks
  }

#pragma omp barrier

#pragma omp for schedule(static) private(i,j)
  for(i=1;i<=n;i++) {                // Loop over vertices
    double heaviest = -1.0;
    partner = 0;
    for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
      int y = edges[j];            // y is the current neighbor of i
//      if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y

        if ((weight[j] < heaviest) || (weight[j] < ws[y]))   // If we already have a better candidate or if y has a better candidate
          continue;					     // Skip to next

        if ((weight[j] == heaviest) && (y < partner))        // Carefull comparison when weights are equal
          continue;

        if ((weight[j] == ws[y]) && (i < s[y]))
          continue;

        heaviest = weight[j];      // Store the best so far
        partner = y;
      
    } // loop over neighbors

    if (heaviest > 0) {
      omp_set_lock(&nlocks[partner]);    // Locking partner
      if (heaviest >= ws[partner]) {
//      if ((heaviest > ws[partner]) || ((heaviest == ws[partner]) && (i>s[partner]))) {
//      If this is true then i will become the new suitor of "partner"
        if (s[partner] != 0) {
     // Add s[partner] (old partner) to local list
     
          left[start+num_left] = s[partner];
          num_left++;
        }

        s[partner] = i;          // i is now the current suitor of s
        ws[partner] = heaviest;  // The weight of the edge (i,partner) 
      }
      else {   // can no longer use the result for node i, must find new partner
      // Add i to the local list
        left[start+num_left] = i;
        num_left++;
      }
      omp_unset_lock(&nlocks[partner]); // Unlocking partner
    }
  } // loop over vertices

  // printf("%d, Remaining vertices %d \n",my_id,num_left);

// Now do the same as above but using the local list


  while(num_left>0) {              // Treat remaining vertices as a stack
    double heaviest = -1.0;
    partner = 0;
    num_left--;                    // Remove the topmost vertex
    i = left[start+num_left];

    for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
      int y = edges[j];            // y is the current neighbor of i
//      if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y
      if ((weight[j] < heaviest) || (weight[j] < ws[y]))
        continue;

      if ((weight[j] == heaviest) && (y < partner))
        continue;

      if ((weight[j] == ws[y]) && (i < s[y]))
        continue;
      heaviest = weight[j];      // Store the best so far
      partner = y;
    } // loop over neighbors

    if (heaviest > 0) {                  // Check if there is a potential partner
      omp_set_lock(&nlocks[partner]);    // Locking partner
      if (heaviest >= ws[partner]) {
//      if ((heaviest > ws[partner]) || ((heaviest == ws[partner]) && (i>s[partner]))) {
        if (s[partner] != 0) {
     // Add s[partner] to local list
          left[start+num_left] = s[partner];
          num_left++;
        }
        s[partner] = i;          // i is now the current suitor of s
        ws[partner] = heaviest;  // The weight of the edge (i,partner) 
      }
      else {   // can no longer use the result for node i, must find new partner
      // Add i to local list
        left[start+num_left] = i;
        num_left++;
      }
      omp_unset_lock(&nlocks[partner]); // Unlocking partner
    }  
  } // while-loop over vertices

#pragma omp barrier
}
