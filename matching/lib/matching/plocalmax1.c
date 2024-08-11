// Round based local max algorithm.
//
// Based on the master thesis of Marcel Birn
// Repeatedly go through edge list and find dominating edges.
// Parallel version
// This time, only using local edge lists. This saves copying back of remaining edges to global list
// but might lead to worse load balance.

void plocalmax1(int n,int m, wed *edgel,int *p,double *ws,int *init,omp_lock_t *nlocks,int *new_m,wed *new_list) {

  int i,j;
  int x,y;
  double w;
  wed *temp;
  // wed *new_list;
  int threads = omp_get_num_threads();
  int my_id = omp_get_thread_num();


  // new_list = (wed *) malloc(sizeof(wed) * m);  // Allocate local memory to store edge list

// Initialize

#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {
    p[i] = 0;          // Set that vertex i is unmatched
    ws[i] = -1.0;      // Set the weight of the best candidate for vertex i to -1.
    omp_init_lock(&nlocks[i]); // Initialize locks
  }

  int tm = 0;  // This is the number of edges for a particular thread

// Copy edges for each thread into local data structure

#pragma omp for schedule(static) private(i)
  for(i=0;i<m;i++) {
    new_list[tm].x = edgel[i].x;
    new_list[tm].y = edgel[i].y;
    new_list[tm].w = edgel[i].w;
    tm++;
  }

#pragma omp barrier

  int r=0;		// Counting the number of rounds
  int old_m = m;	// Used to keep track of how many edges are removed in each round
  int total = m;

  while (total > 0) {	// While there are edges left
/*
    r++;
    if (my_id == threads-1)    // The last thread prints
      printf("Round number %d, remaining edges %d, removed %d \n",r,total,old_m - total);
*/
// First find the best candidate for each vertex

    for(i=0;i<tm;i++) {
      x = new_list[i].x;
      y = new_list[i].y;
      w = new_list[i].w;

      if (w >= ws[x]) {            // Check if x can get a better value with this edge
        omp_set_lock(&nlocks[x]);  // Test must break equality consistently
        if ((w > ws[x]) || ((w == ws[x]) && (y > p[x]))) { // Confirm that x gets a better value with this edge
          ws[x] = w;
          p[x] = y;
        }
        omp_unset_lock(&nlocks[x]);
      }

      if (w >= ws[y])  {           // Check if y can get a better value with this edge
        omp_set_lock(&nlocks[y]);
        if ((w > ws[y]) || ((w == ws[y]) && (x > p[y]))) { // Confirm that y gets a better value with this edge
          ws[y] = w;
          p[y] = x;
        }
        omp_unset_lock(&nlocks[y]);
      }
    }

#pragma omp barrier

// Each vertex is now pointing to its best neighbor
// Go through the edges again and keep those that are incident on two unmatched vertices

    int new_tm = 0;  // This is the number of remaining edges for a particular thread

    for(i=0;i<tm;i++) {
      x = new_list[i].x;
      y = new_list[i].y;

      if ((p[p[x]] != x) && (p[p[y]] != y)) {  // Check if x and y are unmatched

        new_list[new_tm].x = x;            // Store this edge
        new_list[new_tm].y = y;
        new_list[new_tm].w = new_list[i].w;
        new_tm++;

        ws[x] = -1;			// Set that x and y are unmatched
        ws[y] = -1;                     // This might be done by several threads simultaneously
        p[x] = 0;
        p[y] = 0;
      } // end if

    } // end for

//    printf("Thread %d removed %d edges \n",my_id,tm-new_tm);

    new_m[my_id] = new_tm;                  // Store number of remaining edges for each thread
    tm = new_tm;

#pragma omp barrier

// Check if there are edges left to process
    old_m = total;
    total = 0;
    for(i=0;i<threads;i++)
      total += new_m[i];

    if (old_m - total <= 0) {
      printf("Error, removed %d edges \n",old_m-total);
      return;
    }
  } // end while


// Clean up so that check-matching  works correctly

#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {  // Set all unmatched vertices to point to 0
    if (p[p[i]] != i)
      p[i] = 0;
  }
}
