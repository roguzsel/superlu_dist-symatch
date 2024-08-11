// Parallel greedy matching
// Based on traversing the edges followed by a verification step

void edgec1(int n,int m,edge *e,int *p,int *left) {

  int i,j;
  int x,y;
  int nr_m= 0;

  int my_id = omp_get_thread_num();
  int threads = omp_get_num_threads();

  int start = my_id*n/threads;   // Calculate where to start storing elements in "left"

// The start index is not working correctly. Therefore the size of left has been set to 2*n
 
// Start of matching algorithm

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++)   // Set that every node is unmatched
    p[i] = 0;         // Note that p[0] = 0

#pragma omp barrier

#pragma omp for schedule(static) private(i)
  for(i=0;i<m;i++) {                // Loop over edges
    x = e[i].x;                     // Pick up endpoints of the i'th edge
    y = e[i].y;

//    if ((x < 1) || (y<1) || (x>n) || (y>n))
//      printf("%d Has %dth edge = %d,%d \n",my_id,i,x,y);
      
    if ((p[x] == 0) && (p[y] == 0))  {              // Check if both endpoints are available
      p[x]++;                                       // Mark that the x'th vertex is matched
      p[y]++;                                       // Mark that the y'th vertex is matched

//      if (start+nr_m > n)
//        printf("%d is indexing out of bounds %d\n",my_id,start+nr_m);

      left[start+nr_m] = i;                         // Store index of edge for verification
      nr_m++;                                       // Keep track of number of edges matched
   }
  } // loop over edges

// Collect unmatched vertices that still stand a chance of being matched

//  if (start+nr_m > n)
//    printf("%d INDEXING OUT OF BOUNDS!!!: %d \n",my_id,start+nr_m);

#pragma omp barrier

  int un_matched = 0;

#pragma omp for schedule(static) private(i)
  for(i=0;i<nr_m;i++) {                 // Loop over matched edges
    x = e[left[start+i]].x;
    y = e[left[start+i]].y;
    if ((p[x] != 1) || (p[y] != 1)) {   // Check if edge is incident on some other matched edge 
      left[start+un_matched] = left[start+i]; // Store index of edge for later processing
      p[x]--;                           // Remove the edge
      p[y]--;                           // Remove the edge
      un_matched++;
    } 
  } // for i

  if (un_matched != 0) {
    printf("After 1st matching %d has %d edges remaining out of %d \n",my_id,un_matched,m/threads);
  }
  return;
/*

#pragma omp barrier

// int itt = 1;
// Perform remaining matchings
  while (!done) {

//    if (my_id == 0)
//      printf("Iteration %d of the matching \n",itt++);

    for(i=0;i<un_matched;i++) {                   // Loop over remaining unmatched vertices
      int k = left[start+i];                      // Pick up the next unmatched vertex k
      if (p[p[k]] != 0) {                         // Check if vertex k has been matched already
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

    int new_left = 0;
#pragma omp barrier
    done = TRUE;

    for(i=0;i<un_matched;i++) {                 // Loop over remaining vertices
      int k = left[start+i];                    // Pick up vertex k for checking
      if (p[p[k]] != k) {                       // Check if vertex k has been matched
        for(j=ver[k];j<ver[k+1];j++) {          // If not, loop over neighbors of vertex k
          int nv = edges[j];                    // nv is the name of the current neighbor of k
          if (p[p[nv]] != nv) {                 // Check if vertex nv is unmatched
            p[k] = 0;                           // Make sure k is ready for matching
            left[start+new_left] = k;           // Store k in list for later processing
            new_left++;                         // Increase number of unmatched vertices
            break;                              // Continue to next vertex
          }
        } // for j
      } 
    }  // for i

//    if (new_left != 0)
//      printf("Finaly, %d has %d vertices remaining out of %d \n",my_id,new_left,n/threads);
    un_matched = new_left;
#pragma omp barrier
    if (un_matched > 0)
      done = FALSE;
#pragma omp barrier
  }
*/
}
