// Parallel greedy matching
// Based on traversing the edges

void locking1(int n,int m,edge *e,int *p,omp_lock_t *nlocks) {

  int i,j;
  int x,y;

  int my_id = omp_get_thread_num();
  int threads = omp_get_num_threads();

// Start of matching algorithm

//  printf("Starting lock based algorithm \n");

#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {   // Set that every node is unmatched
    p[i] = 0;         
    omp_init_lock(&nlocks[i]); // Initialize locks
  }

#pragma omp for schedule(static) private(i)
  for(i=0;i<m;i++) {                // Loop over vertices
    x = e[i].x;
    y = e[i].y;
    if ((p[x] == 0) && (p[y] == 0))  {              // Check if both endpoints are available
      if (x<y) {
        int t = x;
        x = y;
        y = t;
      }
      omp_set_lock(&nlocks[x]); // locking for x (the larger)
      omp_set_lock(&nlocks[y]); // locking for y (the smaller)
      if ((p[x] == 0) && (p[y] == 0))  {              // Check if both endpoints are still available
        p[x]++;                                       // Mark that vertex is matched
        p[y]++;                                       // Mark that vertex is matched
      }
      omp_unset_lock(&nlocks[x]); // locking for x (the larger)
      omp_unset_lock(&nlocks[y]); // locking for y (the smaller)
    }
  } // loop over edges
#pragma omp barrier
  return;
}
