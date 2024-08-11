// Manne-Bisseling algorithm
//
// Using local lists to hold the unprocessed matched edges
// 
// Parallel version

void MBLocal(int n,int *ver,int *edges,int *p,double *ws,omp_lock_t *nlocks,double *weight,int *locallist,int *num_match,int *glob_list,int *matched,double *weight1) {

  int i,j,k;
  int x,y,z;
  double w;
  wed *temp;
  // wed *new_list;
  int threads = omp_get_num_threads();
  int my_id = omp_get_thread_num();

  int *new_list = (int *) weight1; // Reusing memory as integer

  int *l1 = locallist;
  int *l2 = new_list;
  int *l3;
  int round=0;

  // new_list = (wed *) malloc(sizeof(wed) * m);  // Allocate local memory to store edge list

// Initialize
  p[0] = 0;
#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {
    matched[i] = 0;
  }

// Find the heaviest neighbor
//
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {
    double heaviest = -1.0;
    int partner = 0;
    for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
      int y = edges[j];            // y is the current neighbor of i

      if (weight[j] < heaviest)                          // If the weight of the current edge is lighter than the current best or
        continue;                                        // if we cannot offer y something better, move on.

      if ((weight[j] == heaviest) && (y < partner))      // When equal weight, give priority to highest index
        continue;

      heaviest = weight[j];  
      partner = y;
    }
    p[i] = partner;
  }


// Now every vertex is pointing to its heaviest neighbor
#pragma omp barrier

  int num = 0;
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {
    if ((p[p[i]] == i) && (i < p[i])) {
      matched[i] = true;
      matched[p[i]] = true;
      locallist[num] = i;
      num++;
    }
  }


// Each thread now has a list of candidates

// Calculate how many candidates there are in total

    num_match[my_id] = num;

#pragma omp barrier

    int glob_match = 0;
    for(i=0;i<threads;i++)
      glob_match += num_match[i];
    


  // printf("My_id %d, threads %d \n",my_id,threads);
  //round = 0;

  while (glob_match > 0) {
    round++;
    if (my_id == 0)
      printf("Round %d, global edges %d \n",round,glob_match);

    int new_num = 0;   // The number of matched edges found by this thread in the current iteration

// Loop through the old matched edges

    for(i=0;i<num;i++) {
      int x = l1[i];                 // x was matched to p[x] in the previous round
      for(j=ver[x];j<ver[x+1];j++) { // Loop over neighbors of vertex x
        int y = edges[j];            // y is the current neighbor of x
        if (y == p[x])               // Skip the matching partner of x
          continue;
        if (p[y] == x) {             // y wanted x, must find a new candidate for y
          double heaviest = -1.0;
          int partner = 0;
          for(k=ver[y];k<ver[y+1];k++) { // Loop over neighbors of vertex y
            int z = edges[k];            // z is the current neighbor of y
            if (weight[k] < heaviest)    // If the weight of the current edge is lighter than the current best 
              continue;                  // then move on.
            if (matched[z])              // If z is already matched, ignore it
              continue;
            if ((weight[k] == heaviest) && (y < partner))      // When equal weight, give priority to highest index
              continue;

            heaviest = weight[k];    // z is now the prefered partner of y
            partner = z;
          } // end for k
          p[y] = partner;            // Store which vertex y wants to match with
          if (p[partner] == y) {     // check if y and partner are matched
            matched[y] = true;
            matched[partner] = true;
            l2[new_num] = y;         // Store y for processing in the next round
            new_num++;
          }
        } // end if p[x]
      } // end for j
      x = p[x];     // Now repeat everything for p[x]
      for(j=ver[x];j<ver[x+1];j++) { // Loop over neighbors of vertex x
        int y = edges[j];            // y is the current neighbor of x
        if (y == p[x])               // Skip the partner of x
          continue;
        if (p[y] == x) {             // Must find new candidate for y
          double heaviest = -1.0;
          int partner = 0;
          for(k=ver[y];k<ver[y+1];k++) { // Loop over neighbors of vertex y
            int z = edges[k];            // z is the current neighbor of y
            if (weight[k] < heaviest)    // If the weight of the current edge is lighter than the current best 
              continue;                  // move on.
            if (matched[z])              // If z is already matched, ignore it
              continue;
            if ((weight[k] == heaviest) && (y < partner))      // When equal weight, give priority to highest index
              continue;

            heaviest = weight[k];    // z is now the prefered partner of y
            partner = z;
          } // end for k
          p[y] = partner;            // Store which vertex y wants to match with
          if (p[partner] == y) {     // check if y and partner are matched
            matched[y] = true;
            matched[partner] = true;
            l2[new_num] = y;         // Store y for processing in the next round
            new_num++;
          }
        } // end if p[x]
      } // end for j
    } // end for i

// Calculate the number of new matched edges
    num_match[my_id] = new_num;

#pragma omp barrier

    glob_match = 0;
    for(i=0;i<threads;i++)
      glob_match += num_match[i];

    num = new_num;
// Switch locallists
    l3 = l1;
    l1 = l2;
    l2 = l3;

  } // end while

#pragma omp barrier
  return;

}
