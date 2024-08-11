// Manne-Bisseling algorithm
//
// Using one common list to hold the unprocessed matched edges
// 
// Parallel version

void MB1List(int n,int *ver,int *edges,int *p,double *ws,omp_lock_t *nlocks,double *weight,int *locallist,int *num_match,int *glob_list,int *matched,int *when) {

  int i,j,k;
  int x,y,z;
  double w,heaviest;
  wed *temp;
  // wed *new_list;
  int threads = omp_get_num_threads();
  int my_id = omp_get_thread_num();
  int partner;
  int c,xx;
  int round,kk;
  int my_pos,num;
  int glob_match,xk;


  // new_list = (wed *) malloc(sizeof(wed) * m);  // Allocate local memory to store edge list

// Initialize
  p[0] = 0;
#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {
    matched[i] = false;
    when[i] = 0;
  }

// Find the heaviest neighbor
//
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {
    heaviest = -1.0;
    partner = 0;
    for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
      y = edges[j];            // y is the current neighbor of i

      if (weight[j] < heaviest)                          // If the weight of the current edge is lighter than the current best or
        continue;                                        // if we cannot offer y something better, move on.

      if ((weight[j] == heaviest) && (y < partner))      // When equal weight, give priority to highest index
        continue;

      heaviest = weight[j];  
      partner = y;
    }
    p[i] = partner;	// Note that if there is no candidate then p[i] is set to 0
  }


// Now every vertex is pointing to its heaviest neighbor
#pragma omp barrier

  num = 0;
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {
    // printf("%d is pointing to %d \n",i,p[i]);
    if ((p[p[i]] == i) && (i < p[i])) {      // If we have a matched edge, store the smallest endpoint (only once)
      matched[i] = true;
      matched[p[i]] = true;
      locallist[num] = i;
      num++;
    //  printf("%d and %d are matched \n",i,p[i]);
    }
  }


// Each thread now has a list of candidates

// Calculate where each thread should place its local list in the global edge list
// Do a prefix sum computation on num_match 

    num_match[my_id] = num;  // num_match[] is a global data structure

#pragma omp barrier

    my_pos = 0;
    for(i=0;i<my_id;i++)
      my_pos += num_match[i];

    if (my_id == threads-1) {
      num_match[threads] = my_pos + num;   // This sets the global number of matched edges
    }
#pragma omp barrier

  glob_match = num_match[threads]; // Everybody gets the global number of unprocessed matched edges

  

  round = 0;

  while (glob_match > 0) {
    round++;
/*
    if (my_id == 0) 
       printf("Round %d, matched edges %d \n",round,glob_match);
*/

// Pack the remaining edges back into the global edge list


    for(i=0;i<num;i++) {
      glob_list[i+my_pos] = locallist[i];
    }

    num = 0;

#pragma omp barrier
    xk=0;
    kk=0;


// Do a parallel loop through the matched edges
#pragma omp for schedule(static) private(i)
    for(i=0;i<glob_match;i++) {
      x = glob_list[i];          // x has been matched with p[x]
      kk++; 
      for(j=ver[x];j<ver[x+1];j++) { // Loop over neighbors of vertex x
        xk++;
        y = edges[j];            // y is the current neighbor of x
        if (y == p[x])               // Skip the partner of x
          continue;
        if (p[y] == x) {             // Must find new candidate for y
          heaviest = -1.0;
          partner = 0;
          for(k=ver[y];k<ver[y+1];k++) { // Loop over neighbors of vertex y
            z = edges[k];            // z is the current neighbor of y
            if (weight[k] < heaviest)    // If the weight of the current edge is lighter than the current best or
              continue;                  // if we cannot offer y something better, move on.
            if (matched[z] == true)      // If z is already matched, ignore it
              continue;
            if ((weight[k] == heaviest) && (z < partner))      // When equal weight, give priority to highest index
              continue;

            heaviest = weight[k];    // z is now the prefered partner of y
            partner = z;
          } // end for k
          p[y] = partner;
          // printf("a)%d is now pointing to %d \n",y,partner);
          when[y] = round;  // Mark y with the current round
          locallist[num] = y;  // Store y for later checking
          num++;
        } // end if p[x]
      } // end for j
      x = p[x];     // Now repeat everything for p[x]
      for(j=ver[x];j<ver[x+1];j++) { // Loop over neighbors of vertex x
        xk++;
        y = edges[j];            // y is the current neighbor of x
        if (y == p[x])               // Skip the partner of x
          continue;
        if (p[y] == x) {             // Must find new candidate for y
          heaviest = -1.0;
          partner = 0;
          for(k=ver[y];k<ver[y+1];k++) { // Loop over neighbors of vertex y
            z = edges[k];            // z is the current neighbor of y
            if (weight[k] < heaviest)    // If the weight of the current edge is lighter than the current best or
              continue;                  // if we cannot offer y something better, move on.
            if (matched[z] == true)      // If z is already matched, ignore it
              continue;
            if ((weight[k] == heaviest) && (z < partner))      // When equal weight, give priority to highest index
              continue;

            // printf("%d now has candidate %d of weight %f \n",y,z,weight[k]);
            heaviest = weight[k];    // z is now the prefered partner of y
            partner = z;
          } // end for k
          p[y] = partner;
          // printf("b)%d is now pointing to %d \n",y,partner);
          when[y] = round;
          locallist[num] = y;
          num++;
        } // end if p[x]
      } // end for j
    } // end for i

#pragma omp barrier

// Now pick up the new matched edges, but only do it from one side

    c = 0;             // c is the number of new matched edges
    xx = 0;
    for(i=0;i<num;i++) {
      y = locallist[i];        // Check if y is matched
      partner = p[y];          // y is pointing to partner
      if (p[partner] == y) {   // if partner is pointing to y then they are matched
        xx++;
        if ((y < partner) || (when[partner] != round)) {

          matched[y] = true;
          matched[partner] = true;
          locallist[c] = y;    // store the smaller one for the next round
          c++;
        }
      }
    }

// Calculate the number of new matched edges
    num = c;
    num_match[my_id] = num;

#pragma omp barrier

    my_pos = 0;
    for(i=0;i<my_id;i++)
      my_pos += num_match[i];

    if (my_id == threads-1) {
      num_match[threads] = my_pos + num;
    }
#pragma omp barrier
    glob_match = num_match[threads]; // Everybody gets the global number of unprocessed matched edges
/*
    if ((round ==1) && (glob_match != 12538)) {
      if (my_id == 0) {
        printf("Got %d in round %d, should have been 12538, n=%d\n",glob_match,round,n);
        for(i=1;i<n;i++) {
          if (matched[i])
            printf("%d\n",i);
        }
      }
      #pragma omp barrier
      exit(0);
    }
    if (my_id == 0)
      printf("Round %d: Number of matched edges is %d \n",round,glob_match);
*/

  } // end while

#pragma omp barrier
  return;

}
