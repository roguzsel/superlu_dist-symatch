// weighted algorithm
//
// This is a dummy program to check scalability

void weighted3(int n,int m,int *ver,int *edges,int *p,float *ws,omp_lock_t *nlocks,float *weight,int *left) {


  int i,j;
  int partner;     // Prospective candidate to match with   
  int next_vertex;

  int my_id = omp_get_thread_num();
  int threads = omp_get_num_threads();

  int done,x;

  int sum;
// Start of matching algorithm

//  int start = my_id*n/threads;
  int start = my_id*n/threads;   // Calculate where to start storing elements in "left"

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++)   // Set that every node is unmatched
    p[i] = 0;         // Note that p[0] = 0


#pragma omp barrier

#pragma omp for schedule(static) private(i,j)
  for(i=1;i<=n;i++) {                // Loop over vertices
    if (p[p[i]] != i) {              // Check if vertex i has been matched
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int k = edges[j];            // k is the current neighbor of i
        sum = 0;
        if (p[p[k]] != k) {          // Check if vertex k is unmatched
          for(x=0;x<40;x++)
            sum += x;
            p[i] += (k+sum) % 10;                  // Match i and k
          p[k] = i;
//          break;                     // Continue to next vertex
        }
      }
    } // if
  } // loop over vertices

/*
#pragma omp barrier

  int un_matched = 0;
//  printf("Checking %d vertices \n",n);

#pragma omp for schedule(static) private(i)
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
*/

/*
  if (un_matched != 0)
    printf("Unmatched = %d \n",un_matched);
#pragma omp for schedule(static) private(i)
  for(i=1;i<n;i++) {                // Loop over edges
    sum = edges[i];
    for(j=ver[i];j<ver[i+1];j++) {
      sum += (j+i+edges[j]) % 2;
    }
    edges[i]+= sum;
//      int y = edges[j];            // y is the current neighbor of i
//      if (weight[j]>heaviest) {
//        heaviest = weight[j];      // Store the best so far
//        partner = y;
//      }
//    } // loop over neighbors
  } // loop over vertices

*/

}
