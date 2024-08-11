// Sequential greedy matching
//

void verify(int n,int *ver,int *edges,int *p,int *left) {

  int i,j;
  int nr_union = 0;
  int nr_matched = 0;

  int my_id = omp_get_thread_num();
  int threads = omp_get_num_threads();

  int start = my_id*n/threads;   // Calculate where to start storing elements in "left"

// Start of matching algorithm


#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++)   // Set that every node is unmatched
    p[i] = 0;         // Note that p[0] = 0

#pragma omp barrier

#pragma omp for schedule(static) private(i,j)
  for(i=1;i<=n;i++) {                // Loop over vertices
    if (p[p[i]] != i) {              // Check if vertex i has been matched
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int k = edges[j];            // k is the current neighbor of i
        if (p[p[k]] != k) {          // Check if vertex k is unmatched
          p[i] = k;                  // Match i and k
          p[k] = i;
          break;                     // Continue to next vertex
        }
      }
    } // if
  } // loop over vertices

// Collect unmatched vertices that still stand a chance of being matched


#pragma omp barrier

  int un_matched = 0;
//  printf("Checking %d vertices \n",n);

#pragma omp for schedule(static) private(i,j)
  for(i=1;i<=n;i++) {                // Loop over vertices
    if (p[p[i]] != i) {              // Check if vertex i has been matched
      p[i] = 0;                      // Reset the vertex so that it can be matched
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int k = edges[j];            // k is the current neighbor of i
        if (p[p[k]] != k)   {        // Check if vertex k is unmatched
          left[start+un_matched] = i;
          un_matched++;
          break;                     // Exit the j-loop since i has already been entered
        }
      }
    } 
  } // for i

  if (un_matched != 0) {
    done = false; 
  }
  // printf("After 1st matching %d has %d vertices remaining out of %d \n",my_id,un_matched,n/threads);

#pragma omp barrier

  int itt = 1;
// Perform remaining matchings
  while (!done) {

//    if (my_id == 0)
//      printf("Iteration %d of the matching \n",itt++);

    for(i=0;i<un_matched;i++) {                   // Loop over remaining unmatched vertices
      int k = left[start+i];                      // Pick up the next unmatched vertex k
      if (p[p[k]] != k) {                         // Check if vertex k has been matched already
        for(j=ver[k];j<ver[k+1];j++) {            // Loop over neighbors of vertex k
          int nv = edges[j];                      // nv is the name of the current neighbor
          if (p[p[nv]] != nv)  {                  // Check if nv is matched 
            p[k] = nv;                            // If not, match k and nv
            p[nv] = k;
            break;                                // Continue to next vertex
          }
        } // for j
      } // if
    } // for i

// Now check if there are any vertices that might still be matched
//    if (my_id == 0)
//      printf("Checking if there are any vertices left \n");

    int new_left = 0;
#pragma omp barrier
    done = true;
#pragma omp barrier

    for(i=0;i<un_matched;i++) {                 // Loop over remaining vertices
      int k = left[start+i];                    // Pick up vertex k for checking
      if (p[p[k]] != k) {                       // Check if vertex k has been matched
        p[k] = 0;                               // Make sure k is ready for matching
        for(j=ver[k];j<ver[k+1];j++) {          // If not, loop over neighbors of vertex k
          int nv = edges[j];                    // nv is the name of the current neighbor of k
          if (p[p[nv]] != nv) {                 // Check if vertex nv is unmatched
            left[start+new_left] = k;           // Store k in list for later processing
            new_left++;                         // Increase number of unmatched vertices
            break;                              // Continue to next vertex
          }
        } // for j
      } 
    }  // for i

//    if (my_id == 0)
//      printf("Done checking, found %d vertices \n",new_left);

//    if (new_left != 0)
//      printf("Finaly, %d has %d vertices remaining out of %d \n",my_id,new_left,n/threads);
    un_matched = new_left;

//    if (my_id == 0)
//      printf("Unmatched is now %d \n",un_matched);

#pragma omp barrier
    if (un_matched > 0)
      done = false;
#pragma omp barrier
//    if (my_id == 0)
//      printf("done is now %d \n",done);
  }
//  printf("Exiting \n");
}
