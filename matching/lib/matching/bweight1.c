// parallel weighted algorithm
//
// Computes greedy 2-matching and performs dynamic programming
// Dynamic programming is done in a separate routine 
// One routine for the paths, and one routine for the cycles

//void bweight1(int n,int *ver,int *edges,int *s,double *ws,int *s2, double *ws2,double *weight, int *used,int *match) {
void bweight1(int n,int *ver,int *edges,int *s,double *ws,int *s2, double *ws2,double *weight, int *used,int *match,omp_lock_t *nlocks,int *par,int *path1,double *weight1) {

  int i,j,k;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;
  int count,current;
//  int *path1;
  int l1,l2;
  int next,length1,tmp;
//  double *weight1;
  double cum_weight = 0.0;
  double dyn_prog();
  int skip;
  double two_weight = 0.0;


//  path1 = (int *) malloc(n * sizeof(int));
//  weight1 = (double *) malloc(n * sizeof(double));

// Start of matching algorithm

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0

    s2[i] = 0;           // Set that no node is trying to match with i 
    ws2[i] = 0.0;        // Set the current weight of best suitor edge to 0

    used[i] = false;     // Set that node i has not been used in the dynamic programming part
    match[i] = 0;

    omp_init_lock(&nlocks[i]); // Initialize locks
  }
//  count = 0;


#pragma omp barrier

  i = 1;
  current = i;
  

#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {
    for(k=1;k<3;k++) {
      done = false;
      current = i;
      while (!done) {
//        double heaviest = ws2[current];          // Find the heaviest candidate that does not have a better offer
//        partner = s2[current];                  // No point in trying a partner worse than the best suitor
        double heaviest = -1.0;
        for(j=ver[current];j<ver[current+1];j++) { // Loop over neighbors of the current vertex 
          int y = edges[j];            // y is the neighbor of the current vertex

// Check if w(current,y) is the best so far, and if it is a better option for y
// Also check that current is not already a suitor of y

          if ((weight[j] < heaviest) || (ws2[y] > weight[j])) // Must beat previous heaviest and worst suitor
            continue;
          if ((s[y] == current) || (s2[y] == current))  // Cannot already be pointing to this vertex
            continue;
          if ((weight[j] == heaviest) && (y < partner)) // If equal weight, must have higher index than previous heaviest 
            continue;

//          if ((weight[j] == ws2[y]) && (current < s2[y])) // If equal weight, must have higher index than previous suitor 
//            continue;

          if (weight[j] == ws2[y])  {
            skip = false;
            omp_set_lock(&nlocks[y]);   // Locking y 
              if ((weight[j] == ws2[y]) && (current < s2[y]))            // Must have higher index to beat previous suitor when equal weight
                skip = true;
            omp_unset_lock(&nlocks[y]); // Unlocking y
            if (skip)
              continue;
          }


          heaviest = weight[j];      // Store the weight of the heaviest edge found so far
          partner = y;               // Store the name of the associated neighbor
        } // loop over neighbors

        done = true;
        if (heaviest > 0) {            // True if there is a new partner

          omp_set_lock(&nlocks[partner]);    // Locking partner

          next_vertex = s2[partner];
// Check if current can become the first suitor of partner
          if ((ws[partner] < heaviest) || ((ws[partner] == heaviest) && (s[partner] < current))) {
            s2[partner] = s[partner];    // Shift the current first suitor to become the second suitor
            ws2[partner] = ws[partner];  // Also move the weight
            s[partner] = current;        // current is now the first suitor of s
            ws[partner] = heaviest;      // The weight of the edge (current,partner) 
          }
          else { // current becomes the second suitor of partner
            s2[partner] = current;        // current is now the second suitor of s
            ws2[partner] = heaviest;      // The weight of the edge (current,partner) 
           }

          omp_unset_lock(&nlocks[partner]); // Unlocking partner

          if (next_vertex != 0) {       // True if partner already had a suitor
            current = next_vertex;      // Pick up the old suitor and continue
            done = false;
          }
          //printf("current = %d, s1(%d) = %d, s2(%d) = %d \n",current,partner,s[partner],partner,s2[partner]); 
        }
      } // while loop
    } // for k=1:2
  } // loop over vertices

//  return;  // Code below here should produce a 1-matching by doing DP in parallel on the 2-matching, but this is not ready yet. FM.

#pragma omp barrier
  

// Starting dynamic programming
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    if ((s[i] != 0) && (s2[i] !=0)) // This is not an endpoint, skip it
      continue;

    if ((s[i] == 0) && (s2[i] ==0)) { // This is a singleton, skip it
      used[i] = true;
      continue;
    }

    current = i;

    if (s[i] == 0) {
      next = s2[i];
      weight1[0] = ws2[i];
    }
    else {
      next = s[i];
      weight1[0] = ws[i];
    }
    length1 = 0;

    two_weight += weight1[0];

    while (next != 0) { // Move through the path picking up vertices and edge weights

      path1[length1] = current;
      used[current] = true;
      length1++;
      int tmp = next;
      if (s[next] != current) {
        next = s[next];
        weight1[length1] = ws[tmp];
      }
      else {
        next = s2[next];
        weight1[length1] = ws2[tmp];
      }
      two_weight += weight1[length1];
      current = tmp;
    }
    path1[length1] = current; // Add the final vertex
    used[current] = true;

    if (path1[0] < path1[length1])    // Make sure DP only starts from one end of the path
      dyn_prog(length1,path1,weight1,match);
  }


#pragma omp barrier
  
// Everything above should work, Fredrik 12.02.2014 

// Now look for the cycles
  int nr_cycle = 0;
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    path1[0] = i;
    weight1[0] = ws[i];
    tmp = i;
    current = s[i]; // Initially start by moving to the right
    length1 = 0;
    used[i] = true;

    two_weight += weight1[0];

// Now tmp and current point to two consecutive nodes on the cycle, tmp has been processed.
//
// // Pick up rest of cycle
    while (current != i) {

      length1++;
      path1[length1] = current;    // Put current in the list
      used[current] = true;
      if (s[current] == tmp) { // Moving next forward to v(i+1)
        next = s2[current];
        weight1[length1] = ws2[current];  // Picking up weight of edge (current,next)
      }
      else {
        next = s[current];
        weight1[length1] = ws[current];
      }
      two_weight += weight1[length1];

      tmp = current;
      current = next;

    } // End while

// Now do dynamic programming 
    cycle_dyn_prog(length1,path1,weight1,match);

  } // End of loop over vertices that looks for unprocessed cycles

#pragma omp barrier

  for(i=1;i<=n;i++) {
    if (match[i] == 0) {
      double best = -1.0;
      for(k=ver[i];k<ver[i+1];k++) {
        int y = edges[k];
        if ((match[y] == 0) && (weight[k] > best)) {
          best = weight[k];
          partner = y;
        }
      }
      if (best > 0.0) {
        match[i] = partner;
        match[partner] = i;
      }
    }
  }

  free(path1);
  free(weight1);

  printf("Bmatching: Cost of 2-matching %f \n",two_weight);
}

