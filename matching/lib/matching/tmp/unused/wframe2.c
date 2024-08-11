#include "wframe2.h"	// Include files and define variables

#include "seq1.c" 		// Sequential unweighted greedy matching,  vertex based 
#include "eseq1.c"		// Sequential unweighted greedy matching,  edge based 

#include "wseq1.c" 		// Sequential weighted localy greedy matching, vertex based
#include "wseq2.c" 		// Sequential weighted localy greedy matching, vertex based, but in random order
#include "wseq3.c" 		// Same as above, but code has been optimized with prefetching directives

#include "ver1.c"               // Parallel unweighted matching, vertex based using verification to check correctness
#include "edge1.c"              // Parallel unweighted matching, edge based using verification to check correctness
#include "edge2.c"              // Optimized version of edge1
#include "lock1.c"              // Parallel unweighted matching, edge based using locking to assure correctness
#include "weight1.c"            // First version of parallel suitor code, not very efficient
#include "weight2.c"            // Final version of parallel suitor code

#include "sweight.c" 		// Suitor based weighted matching. Using a for-loop as outer loop
#include "sweight1.c" 		// Suitor based weighted matching. Using sorted adjacency lists
#include "sweight7.c" 		// Suitor based weighted matching. Optimized loop.
#include "sweight3.c" 		// Suitor based weighted matching. Using a while-loop as outer loop
#include "sweight4.c" 		// Two rounds of suitor based matching, followed by dyn. prog 
#include "sweight5.c" 		// Same as above, but now with separate routine for dyn. prog
#include "sweight6.c" 		// Two rounds of localy greedy matching (wseq), followed by dyn. prog
#include "sweight8.c" 		// Same as sweight5, but now with separate routine for cyclic dyn. prog
#include "sweight9.c" 		// Same as sweight6, but now with separate routine for cyclic dyn. prog
#include "sweight10.c" 		// Two level algorithm with DP, assuming sorted neightbor lists
#include "sweight11.c" 		// Suitor based algorithm with compression of neighbor lists 
#include "pdp.c" 		// Two level parallel suitor algorithm with DP
#include "spdp.c" 		// Two level parallel suitor algorithm with DP
#include "psweight1.c" 		// Parallel suitor based weighted matching. Using sorted adjacency lists

#include "gpa.c" 		// Create paths and cycles by adding in edges (by decreasing weight), followed by dyn. prog
#include "pth1.c" 		// 2 rounds of Path growing algorithm followed by dynamic programming
#include "path.c" 		// Path growing algorithm followed by dynamic programming

#include "wsort.c"		// Sorts edge lists of a graph by decreasing weight. Uses insertion sort
#include "mmio.c"               // Matrix market routines for reading files

// Routine for comparing the weight of two edges
// Used for sorting edges by increasing weight in qsort()

static int compedge(const void *m1, const void *m2) {
    wed *v1 = (wed *) m1; 
    wed *v2 = (wed *) m2; 
    return (int) (v1->w < v2->w);
}

static int compneig(const void *m1, const void *m2) {
    neig *v1 = (neig *) m1; 
    neig *v2 = (neig *) m2; 
    return (int) (v1->w < v2->w);
}

int main(int argc, char *argv[]) {

// Primary graph data structure, using compressed neighbor lists.
// ver stores pointers into the edges list. Note that numbering of both vertices and edges both start at 1

  int *ver;          // Pointers to edge lists
  int *edges;        // Edge lists
  double *weight;    // Corresponding edge weights

// Edge lists where each neighbor list has been sorted by decreasing weight

  int *s_edges;      // Sorted edge lists
  double *s_weight;  // Sorted corresponding edge weights

// Vertex lists used in various algorithms

  int *next;         
  int *used;        // Keeps track of which vertices has been used in the dynamic programming

// Vertex lists used for holding the matching

  int *p;	     // If p(i) = j then i is matched with j, this is used in all the algorithms
  int *p1;           // Same as p, but used for first round of two round algorithms
  int *p2;           // Same as p, but used for second round of two round algorithms

// Vertex lists for suitor algorithm

  double *ws1;       // If ws1(i) = w, then the weight of the best suitor of i is w
  double *ws2;       // Same as ws1, but used in two round algorithm

// Raw edge lists in the order it was read from file

  edge *e;           // List of all edges, no weight, one struct for each edge
  wed *we;	     // Same as e, but now with weights
  neig *swe;	     // List of neighbor + weight 

  int read_graph(),get_input(),allocate_memory();
  double cost_matching();
  void store_value();

  int i,j;          // Loop indices
  double mt1,mt2;   // Used for timing of each algorithm


  if (!get_input(argc,argv,&n_graphs,&n_runs,&n_conf,&conf,&name))	// Get input and check that it is in order
    return(false);

  if (!prepare_output(&wf,n_conf,conf)) {				// Prepare output file
    return(false);
  }

  if (!allocate_memory(n_graphs,&n,&m,&timer,&cost,&p_timer,&p_cost,n_conf)) {		// Allocating memory
    return(false);
  }

  printf("Starting main loop \n");

// Run the main loop, processing each graph, first sequentially and then in parallel

  for(i=0;i<n_graphs;i++) {

    printf("\n");
    printf("***********************************************************\n");
    printf("Starting to read graph number %d: %s \n",i,name[i]);
    printf("***********************************************************\n");
    printf("\n");

// Reading the i'th graph

    if (!read_graph(&(n[i]),&(m[i]),name[i],&ver,&edges,&e,&weight,&we,&swe)) {
      printf("Problem reading graph %s \n",name[i]);
      return(false);
    }

    int xxx = ver[n[i]+1];
    printf("Asking for %d integers \n",xxx);
    s_edges = (int *) malloc(sizeof(int)*xxx);
    printf("Got %d integers \n",xxx);
    if (s_edges == NULL) {
      printf("Unable to allocate memory for s_edges in wframe2 \n");
      return;
    }
    printf("Asking for %d doubles \n",ver[n[i]+1]);
    s_weight = (double *) malloc(sizeof(double)*(ver[n[i]+1]));
    if (s_weight == NULL) {
      printf("Unable to allocate memory for s_weight in wframe2 \n");
      return;
    }
    printf("Got them \n");
    


    if (!allocate_graph_memory(n[i],&next,&used,&p,&p1,&p2)) {			// Allocating memory specifically for this graph
      return(false);
    }


// List of best local suitor for each vertex in the weighted algorithm
// Why is this not in allocate_graph_memory?

    ws1 = (double *) malloc(sizeof(double)*((*n)+1));
    ws2 = (double *) malloc(sizeof(double)*((*n)+1));
    if ((ws1 == NULL) || (ws2 == NULL)) {
      printf("Unable to allocate space for ws-arrays \n");
      return(false);
    }


// Run sequential code 

    printf("\n*******************************\n");
    printf("**** Sequential algorithms ****\n");
    printf("*******************************\n");


    printf("\n**** Unweighted algorithms ****\n");
    int l;
    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      seq1(n[i],ver,edges,p);	// Sequential unweighted greedy vertex based matching
      mt2 = omp_get_wtime(); 
      if ((timer[0][i] < 0) || (mt2-mt1 < timer[0][i]))
        timer[0][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("seq1() did not produce a legal matching \n");
      exit(false);
    }

    printf("Sequential unweighted vertex based greedy matching took %lf seconds\n",timer[0][i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      eseq1(n[i],m[i],e,p);	// Sequential unweighted greedy edge based matching
      mt2 = omp_get_wtime(); 
      if ((timer[1][i] < 0) || (mt2-mt1 < timer[1][i]))
        timer[1][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("eseq1() did not produce a legal matching \n");
      exit(false);
    }

    printf("Sequential unweighted edge based greedy matching took %lf seconds\n",timer[1][i]);
    printf("\n");

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      wseq1(n[i],ver,edges,p,weight); // Sequential weighted localy greedy matching, vertex based
      mt2 = omp_get_wtime(); 
      if ((timer[2][i] < 0) || (mt2-mt1 < timer[2][i]))
        timer[2][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("wseq1() did not produce a legal matching \n");
      exit(false);
    }

    cost[2][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("\n**** Weighted localy greedy algorithms ****\n");
    printf("Sequential local maximal weight matching took %lf seconds, cost=%e\n",timer[2][i],cost[2][i]);

// Compute a random ordering of the vertices
// This is used to traverse the vertices of the graph in a random order

    if (!random_order(n[i],next)) {
      printf("Problem in random_order() \n");
      exit(false);
    }

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      wseq2(n[i],ver,edges,p,weight,next); // Sequential weighted localy greedy matching, vertex based, but in random order
      mt2 = omp_get_wtime(); 
      if ((timer[3][i] < 0) || (mt2-mt1 < timer[3][i]))
        timer[3][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("wseq2() did not produce a legal matching \n");
      exit(false);
    }

    cost[3][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching


// Do the same as above but now in-order

    for(l=1;l<n[i];l++)
      next[l] = l;

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime();
      wseq2(n[i],ver,edges,p,weight,next); // Sequential weighted localy greedy matching, vertex based
      mt2 = omp_get_wtime();
      if ((timer[21][i] < 0) || (mt2-mt1 < timer[21][i]))
        timer[21][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("wseq2() did not produce a legal matching \n");
      exit(false);
    }

    cost[21][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching


    printf("Random order local maximal weight matching took %lf seconds, cost=%e\n",timer[3][i],cost[3][i]);

/*
    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      wseq3(n[i],ver,edges,p,weight,next); // Same as above, but code has been optimized with prefetching directives
      mt2 = omp_get_wtime(); 
      if ((timer[4][i] < 0) || (mt2-mt1 < timer[4][i]))
        timer[4][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("wseq3() did not produce a legal matching \n");
      exit(false);
    }
    cost[4][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Random order local maximal weight matching ** With Prefetching ** took %lf seconds, cost=%e\n",timer[4][i],cost[4][i]);
    printf("\n");
*/

    printf("\n**** Suitor based algorithms ****\n");
    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight(n[i],ver,edges,p,ws1,weight,p1);     // Suitor based weighted matching. Using a for-loop as outer loop
      mt2 = omp_get_wtime(); 
      if ((timer[5][i] < 0) || (mt2-mt1 < timer[5][i]))
        timer[5][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight() did not produce a legal matching \n");
      exit(false);
    }

    cost[5][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Sequential suitor-based weighted matching took %lf seconds, cost=%e\n",timer[5][i],cost[5][i]);

/*

    int *nn = (int *) malloc((n[i]+3)*sizeof(int));

    for(l=0;l<n_runs;l++) {

      for(j=0;j<=n[i]+1;j++) 
        nn[j] = ver[j];

      mt1 = omp_get_wtime(); 
      sweight11(n[i],ver,edges,p,ws1,weight,nn);     // Suitor based weighted matching. Using a for-loop as outer loop
      mt2 = omp_get_wtime(); 
      if ((timer[19][i] < 0) || (mt2-mt1 < timer[19][i]))
        timer[19][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight11() did not produce a legal matching \n");
      exit(false);
    }

    cost = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Sequential suitor-based weighted matching with compressed neighbor lists took %lf seconds, cost=%e\n",timer[19][i],cost);
*/


    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight3(n[i],ver,edges,p,ws1,weight);   // Suitor based weighted matching. Using a while-loop as outer loop
      mt2 = omp_get_wtime(); 
      if ((timer[6][i] < 0) || (mt2-mt1 < timer[6][i]))
        timer[6][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight3() did not produce a legal matching \n");
      exit(false);
    }

    cost[6][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Optimized suitor-based sequential weighted matching took %lf seconds, cost=%e\n",timer[6][i],cost[6][i]);


// Copy the neighbor list + weights into list of structs (needed for sorting using qsort)

    for(j=0;j<ver[n[i]+1];j++) {
      swe[j].w = weight[j];
      swe[j].x = edges[j];
    }

    printf("Starting sorting \n");

    mt1 = omp_get_wtime(); 
//    wsort(n[i],ver,edges,weight,&s_edges,&s_weight);  // Sort the edge list of each vertex by decreasing weight
    for(j=1;j<=n[i];j++) {
      qsort(&(swe[ver[j]]),ver[j+1]-ver[j],sizeof(neig),compneig); 		 // Sort the edge list of each vertex by decreasing weight
    }
    mt2 = omp_get_wtime(); 
    timer[20][i] = mt2-mt1;

    printf("Edge list sorting took %lf seconds\n",timer[20][i]);

/*    int xxx = ver[n[i]+1];
    printf("Asking for %d integers \n",xxx);
    s_edges = (int *) malloc(sizeof(int)*xxx);
    printf("Got %d integers \n",xxx);
*/

//    printf("Asking for %d integers \n",ver[n[i]+1]);
 
//    s_edges = (int *) malloc(sizeof(int)*xxx);
    printf("Got them \n");
/*
    if (s_edges == NULL) {
      printf("Unable to allocate memory for s_edges in wframe2 \n");
      return;
    }
    printf("Asking for %d doubles \n",ver[n[i]+1]);
    s_weight = (double *) malloc(sizeof(double)*(ver[n[i]+1]));
    if (s_weight == NULL) {
      printf("Unable to allocate memory for s_weight in wframe2 \n");
      return;
    }

    printf("Got them \n");
*/
    printf("Filling up s_edges and s_weight \n");

    for(j=0;j<ver[n[i]+1];j++) {
      s_edges[j] = swe[j].x;
      s_weight[j] = swe[j].w;
    }


    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight1(n[i],ver,s_edges,p,ws1,s_weight,next);
      mt2 = omp_get_wtime(); 
      if ((timer[13][i] < 0) || (mt2-mt1 < timer[13][i]))
        timer[13][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight1() did not produce a legal matching \n");
      exit(false);
    }
    cost[13][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Sequential suitor-based weighted matching, using sorted edge lists took %lf seconds, cost=%e\n",timer[13][i],cost[13][i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight7(n[i],ver,s_edges,p,ws1,s_weight,next);  // Suitor based code, using sorted adjacency lists
      mt2 = omp_get_wtime(); 
      if ((timer[14][i] < 0) || (mt2-mt1 < timer[14][i]))
        timer[14][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight7() did not produce a legal matching \n");
      exit(false);
    }

    cost[14][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Optimized sequential suitor-based weighted matching, using sorted edge lists took %lf seconds, cost= %e\n",timer[14][i],cost[14][i]);
    printf("\n");


// Next follows the algorithms using dynamic programming for post-processing

    printf("\n**** Algorithms using DP ****\n");
/*
    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight4(n[i],ver,edges,p,ws1,p2,ws2,weight,used);      // Two rounds of suitor based matching, followed by dyn. prog 
      mt2 = omp_get_wtime(); 				     // This does not produce a complete matching, only the weight
      if ((timer[7][i] < 0) || (mt2-mt1 < timer[7][i]))
        timer[7][i] = mt2-mt1;
    }

    printf("Two round sequential weighted matching with dynamic programming took %lf seconds\n",timer[7][i]);
*/

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight5(n[i],ver,edges,p1,ws1,p2,ws2,weight,used,p);      // Same as above, but now with separate routine for dyn. prog
      mt2 = omp_get_wtime(); 
      if ((timer[8][i] < 0) || (mt2-mt1 < timer[8][i]))
        timer[8][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight5() did not produce a legal matching \n");
      exit(false);
    }

    cost[8][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Two round sequential weighted matching with separate dynamic programming took %lf seconds, cost= %e \n",timer[8][i],cost[8][i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight8(n[i],ver,edges,p1,ws1,p2,ws2,weight,used,p);      // Same as above, but now with separate routine for cycle dyn. prog
      mt2 = omp_get_wtime(); 
      if ((timer[15][i] < 0) || (mt2-mt1 < timer[15][i]))
        timer[15][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight8() did not produce a legal matching \n");
      exit(false);
    }

    cost[15][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Two round sequential weighted matching with separate cycle dynamic programming took %lf seconds, cost= %e \n",timer[15][i],cost[15][i]);


    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight6(n[i],ver,edges,p1,ws1,p2,ws2,weight,used,p);      // Two rounds of localy greedy matching (wseq), followed by dyn. prog
      mt2 = omp_get_wtime(); 
      if ((timer[9][i] < 0) || (mt2-mt1 < timer[9][i]))
        timer[9][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight6() did not produce a legal matching \n");
      exit(false);
    }

    cost[9][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Two round greedy weighted matching with separate dynamic programming took %lf seconds, cost=%e \n",timer[9][i],cost[9][i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime();
      sweight9(n[i],ver,edges,p1,ws1,p2,ws2,weight,used,p);      // Two rounds of localy greedy matching (wseq), followed by dyn. prog
      mt2 = omp_get_wtime();
      if ((timer[16][i] < 0) || (mt2-mt1 < timer[16][i]))
        timer[16][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight9() did not produce a legal matching \n");
      exit(false);
    }

    cost[16][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Two round greedy weighted matching with separate cyclic dynamic programming took %lf seconds, cost=%e \n",timer[16][i],cost[16][i]);


    mt1 = omp_get_wtime(); 
    qsort(we,m[i],sizeof(wed),compedge);	// Sort the edges by decreasing weight
    mt2 = omp_get_wtime(); 
    timer[10][i] = mt2 - mt1;


    for(l=0;l<m[i];l++) {  			// Check that the edges are sorted
      if (we[l].w < we[l+1].w) {
        printf("Data is not sorted in descending order \n");
        printf("%d =  %lf   < %d =  %lf \n",l,we[l].w,l+1,we[l+1].w);
        exit(0);
      }
    }


    printf("Time to sort the edges is %lf seconds \n",timer[10][i]);


    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      gpa(n[i],m[i],we,ver,edges,weight,p);	// Create paths and cycles by adding in edges (by decreasing weight), followed by dyn. prog
      mt2 = omp_get_wtime(); 
      if ((timer[11][i] < 0) || (mt2-mt1 < timer[11][i]))
        timer[11][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("GPA() did not produce a legal matching \n");
      exit(false);
    }

    cost[11][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("GPA algorithm with dynamic programming took %lf seconds, cost=%e\n",timer[11][i],cost[11][i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      path(n[i],ver,edges,weight,p);		// Path growing algorithm followed by dynamic programming
      mt2 = omp_get_wtime(); 
      if ((timer[12][i] < 0) || (mt2-mt1 < timer[12][i]))
        timer[12][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("path() did not produce a legal matching \n");
      exit(false);
    }

    cost[12][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Path growing algorithm with dynamic programming took %lf seconds, cost= %e\n",timer[12][i],cost[12][i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      pth1(n[i],ver,edges,weight,p);		// Path growing algorithm followed by dynamic programming
      mt2 = omp_get_wtime(); 
      if ((timer[17][i] < 0) || (mt2-mt1 < timer[17][i]))
        timer[17][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("pth1() did not produce a legal matching \n");
      exit(false);
    }

    cost[17][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("2-level Path growing algorithm followed by dynamic programming took %lf seconds, cost= %e\n",timer[17][i],cost[17][i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime();
      sweight10(n[i],ver,s_edges,p1,ws1,p2,ws2,s_weight,next,used,p);
      mt2 = omp_get_wtime();
      if ((timer[18][i] < 0) || (mt2-mt1 < timer[18][i]))
        timer[18][i] = mt2-mt1;
    }

    if (!verify_matching(n[i],ver,edges,p)) {
      printf("sweight10() did not produce a legal matching \n");
      exit(false);
    }
    cost[18][i] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching

    printf("Sequential 2 round suitor-based weighted matching + DP, using sorted edge lists took %lf seconds, cost=%e\n",timer[18][i],cost[18][i]);


// Allocate space for locks

    omp_lock_t *nlocks;
    nlocks = (omp_lock_t *) malloc((n[i]+1)*sizeof(omp_lock_t));

    printf("\n*****************************\n");
    printf("**** Parallel algorithms ****\n");
    printf("*****************************\n");

// Run parallel algorithm
// One iteration for each thread configuration
    for(j=0;j<n_conf;j++) {
      printf("\nConfiguration %d, using %d threads \n",j,conf[j]);
      omp_set_num_threads(conf[j]);  // Set number of threads in this configuration

#pragma omp parallel 
      {
        int k;
        int threads = omp_get_num_threads();
        int my_id = omp_get_thread_num();
      

        if (threads != conf[j])
          printf("***** Did not get the correct numer of threads! Wanted %d, got %d *****\n",conf[i],threads);

// Now we are ready to run the parallel algorithms
// Run as many times as required, only store timings for best run

#pragma omp barrier

// Time how long time it takes to sort the individual edge lists

// Copy the neighbor list + weights into list of structs (needed for sorting using qsort)
        int b;
#pragma omp for schedule(static) private(b)
        for(b=0;b<ver[n[i]+1];b++) {
          swe[b].w = weight[b];
          swe[b].x = edges[b];
        }

#pragma omp barrier
#pragma omp master
        { mt1 = omp_get_wtime(); }

#pragma omp for schedule(static) private(b)
        for(b=1;b<=n[i];b++) {
          qsort(&(swe[ver[b]]),ver[b+1]-ver[b],sizeof(neig),compneig);               // Sort the edge list of each vertex by decreasing weight
        }

#pragma omp barrier
#pragma omp master
        { mt2 = omp_get_wtime();
          printf("Edge list sorting took %lf seconds\n",mt2-mt1);
        }



#pragma omp barrier

        for(k=0;k<n_runs;k++) {

// First, unweighted greedy matching algorithms

#pragma omp barrier
#pragma omp master
          { // printf("verify\n");
            mt1 = omp_get_wtime(); }

          verify(n[i],ver,edges,p,next);   // Run vertex based verification code

#pragma omp master
          { mt2 = omp_get_wtime();
            if ((p_timer[0][i][j] < 0) || (mt2-mt1 < p_timer[0][i][j]))
              p_timer[0][i][j] = mt2-mt1;

            if (!verify_matching(n[i],ver,edges,p)) {
              printf("verify() did not produce a legal matching \n");
              exit(false);
            }
          }

#pragma omp barrier



#pragma omp master
          { // printf("edgec1\n");
            mt1 = omp_get_wtime(); }

          edgec1(n[i],m[i],e,p,next);   // Run edge based verification code

#pragma omp master
          { mt2 = omp_get_wtime(); 
            if ((p_timer[1][i][j] < 0) || (mt2-mt1 < p_timer[1][i][j]))
              p_timer[1][i][j] = mt2-mt1;

// Note that the current version of edgec does not return the actual matching

          }
#pragma omp barrier

// Running optimized edge based verification code
/*
#pragma omp master
          {  mt1 = omp_get_wtime(); }
             edgec2(n[i],m[i],e,p,next);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t1_edge[i][j] < 0) || (mt2-mt1 < t1_edge[i][j]))
               t1_edge[i][j] = mt2-mt1;
          }

#pragma omp barrier

*/



#pragma omp master
          { // printf("locking1\n");
            mt1 = omp_get_wtime(); }

          locking1(n[i],m[i],e,p,nlocks);  // Run edge based locking code

#pragma omp master
          { mt2 = omp_get_wtime(); 
            if ((p_timer[2][i][j] < 0) || (mt2-mt1 < p_timer[2][i][j]))
              p_timer[2][i][j] = mt2-mt1;

// Note the current version of locking1 does not return the actual matching
          }

#pragma omp barrier



// Running weighted algorithms
/*
#pragma omp master
          { mt1 = omp_get_wtime(); }
          weighted(n[i],ver,edges,p,ws1,next,nlocks,weight);

#pragma omp master
          { mt2 = omp_get_wtime(); 
            if ((p_timer[3][i][j] < 0) || (mt2-mt1 < p_timer[3][i][j]))
              p_timer[3][i][j] = mt2-mt1;

            if (!verify_matching(n[i],ver,edges,p)) {
              printf("weighted() did not produce a legal matching \n");
              exit(false);
            }
          }
*/

// Running second weighted algorithm

#pragma omp barrier

#pragma omp master
          { // printf("weighted2\n");
            mt1 = omp_get_wtime(); }

          weighted2(n[i],ver,edges,p,ws1,nlocks,weight);

#pragma omp master
          { mt2 = omp_get_wtime(); 
            if ((p_timer[4][i][j] < 0) || (mt2-mt1 < p_timer[4][i][j]))
              p_timer[4][i][j] = mt2-mt1;

            if (!verify_matching(n[i],ver,edges,p)) {
              printf("weighted2() did not produce a legal matching \n");
              exit(false);
            }

            p_cost[4][i][j] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching
          }

#pragma omp barrier

// Running parallel suitor algorithm with sorted adjacency lists

#pragma omp barrier

#pragma omp master
          { // printf("psweight1\n");
            mt1 = omp_get_wtime(); }

          psweight1(n[i],ver,s_edges,p,ws1,s_weight,next,nlocks);

#pragma omp master
          { mt2 = omp_get_wtime();
            if ((p_timer[6][i][j] < 0) || (mt2-mt1 < p_timer[6][i][j]))
              p_timer[6][i][j] = mt2-mt1;

            if (!verify_matching(n[i],ver,edges,p)) {
              printf("psweight1() did not produce a legal matching \n");
              exit(false);
            }

            p_cost[6][i][j] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching
          }

#pragma omp barrier


// Running 2 stage weighted algorithm with dynamic programming

#pragma omp barrier

#pragma omp master
          { // printf("pdp\n");
            mt1 = omp_get_wtime(); }

               // pdp(n[i],ver,edges,p,ws1,nlocks,weight);
          pdp(n[i],ver,edges,p1,ws1,p2,ws2,weight,used,p,nlocks,next);  


#pragma omp master
          { mt2 = omp_get_wtime(); 
            if ((p_timer[5][i][j] < 0) || (mt2-mt1 < p_timer[5][i][j]))
              p_timer[5][i][j] = mt2-mt1;

            if (!verify_matching(n[i],ver,edges,p)) {
              printf("pdp() did not produce a legal matching \n");
              exit(false);
            }

            p_cost[5][i][j] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching
          }

#pragma omp barrier

// Running 2 stage weighted algorithm assuming sorted edge lists followed by dynamic programming

#pragma omp barrier

#pragma omp master
          { // printf("spdp\n");
            mt1 = omp_get_wtime(); }

               // pdp(n[i],ver,edges,p,ws1,nlocks,weight);
          spdp(n[i],ver,s_edges,p1,ws1,p2,ws2,s_weight,used,p,nlocks,next);


#pragma omp master
          { mt2 = omp_get_wtime();
            if ((p_timer[7][i][j] < 0) || (mt2-mt1 < p_timer[7][i][j]))
              p_timer[7][i][j] = mt2-mt1;

            if (!verify_matching(n[i],ver,edges,p)) {
              printf("spdp() did not produce a legal matching \n");
              exit(false);
            }

            p_cost[7][i][j] = cost_matching(n[i],ver,edges,weight,p);  // Get the cost of the matching
          }

#pragma omp barrier


        } // End of k iterations over same configuration 


       
#pragma omp master
        { 
          printf("\n**** Unweighted algorithms ****\n");
          printf("Vertex code took %lf seconds, absolute speedup = %4.2lf, relative speedup = %4.2lf \n",p_timer[0][i][j],
                  timer[0][i]/p_timer[0][i][j],p_timer[0][i][0]/p_timer[0][i][j]);
          printf("Edge code took %lf seconds, absolute speedup = %4.2lf, relative speedup = %4.2lf \n",p_timer[1][i][j],
                  timer[1][i]/p_timer[1][i][j],p_timer[1][i][0]/p_timer[1][i][j]);
          // printf("** Optimized edge code using %d threads took %lf seconds, absolute speedup = %4.2lf \n",conf[j],t1_edge[i][j],t_bse[i]/t1_edge[i][j]);
          printf("Locking code took %lf seconds, absolute speedup = %4.2lf, relative speedup %4.2lf \n",p_timer[2][i][j],
                  timer[1][i]/p_timer[2][i][j],p_timer[2][i][0]/p_timer[2][i][j]);

          printf("\n**** Weighted algorithms ****\n");
//         printf("Weighted code took %lf seconds, absolute speedup = %4.2lf, relative speedup = %4.2lf  \n",p_timer[3][i][j],
//                  timer[6][i]/p_timer[3][i][j],p_timer[3][i][0]/p_timer[3][i][j]);
          printf("Parallel suitor algorithm took %lf seconds, absolute speedup = %4.2lf, relative speedup = %4.2lf, cost = %e \n",p_timer[4][i][j],
                  timer[6][i]/p_timer[4][i][j],p_timer[4][i][0]/p_timer[4][i][j],p_cost[4][i][j]);

          printf("Parallel suitor algorithm with sorted lists took %lf seconds, absolute speedup = %4.2lf, relative speedup = %4.2lf, cost = %e \n",
                  p_timer[6][i][j],timer[13][i]/p_timer[6][i][j],p_timer[6][i][0]/p_timer[6][i][j],p_cost[6][i][j]);

          printf("Parallel DP took %lf seconds, absolute speedup = %4.2lf, relative speedup = %4.2lf, cost = %e \n",p_timer[5][i][j],
                  timer[15][i]/p_timer[5][i][j],p_timer[5][i][0]/p_timer[5][i][j],p_cost[5][i][j]);

          printf("Parallel DP with sorted lists took %lf seconds, absolute speedup = %4.2lf, relative speedup = %4.2lf, cost = %e \n",
                  p_timer[7][i][j],timer[18][i]/p_timer[7][i][j],p_timer[7][i][0]/p_timer[7][i][j],p_cost[7][i][j]);
        }

      }  // End of parallel section


    }  // End of different thread configurations

    printf("Freeing up space \n");
    free(next);
    free(nlocks);   // Free up space used for locks
    free(ver);
    free(edges);
    free(e);
    free(weight);
    free(we);
    free(used);
    free(p);
    free(p1);
    free(p2);
    free(ws1);
    free(ws2);
    free(s_edges);
    free(s_weight);
    printf("Done freeing space\n");


  } // Loop over graphs

// Print results to file
//
// name[i] 	Name of file i (including the .mtx extension)



  printf("Printing to file \n");
  char *tc;         // Name of current graph without the .mtx extension
  

  fprintf(wf,"name = char(");       // Print the names of the graphs
  for(i=0;i<n_graphs;i++)  {
    tc = strtok(name[i],"."); 		// Remove .mtx extension from name
    fprintf(wf,"'%s' ",tc);
    if (i != n_graphs-1)
      fprintf(wf,",");
  }
  fprintf(wf,");\n");

  fprintf(wf,"n = [");       // Print the number of vertices of the graphs
  for(i=0;i<n_graphs;i++)  {
    fprintf(wf,"%d ",n[i]);
  }
  fprintf(wf,"];\n");

  fprintf(wf,"nz = [");       // Print the number of non-zeros of the graphs
  for(i=0;i<n_graphs;i++)  {
    fprintf(wf,"%d ",m[i]);
  }
  fprintf(wf,"];\n");

// Sequential timings
// ******************
// timer[0][]	Sequential unweighted greedy vertex based matching
// timer[1][]	Sequential unweighted greedy edge based matching
// timer[3][]	Sequential weighted localy greedy matching, vertex based, but in random order
// timer[4][]	Same as above, but code has been optimized with prefetching directives
// timer[5][] 	Suitor based weighted matching. Using a for-loop as outer loop
// timer[6][]	Suitor based weighted matching. Using a while-loop as outer loop
// timer[7][]	Two rounds of suitor based matching, followed by dyn. prog, only weight, no matching
// timer[8][]	Same as above, but now with separate routine for dyn. prog, will give matching.
// timer[9][]	Two rounds of localy greedy matching (wseq), followed by dyn. prog
// timer[10][]	Sort all of the edges by decreasing weight
// timer[11][]	GPA algorithm with dynamic programming
// timer[12][]	Path growing algorithm followed by dynamic programming
// timer[13][]	Suitor-based weighted matching, using sorted edge lists
// timer[14][] 	Suitor based code, using sorted adjacency lists, "optimized"
// timer[15][]	Two round suitor code with separate routine for cycle dynamic programming
// timer[16][]	Two rounds of localy greedy matching (wseq), followed by separate cyclic dyn. prog, compare with [9]
// timer[17][]	2-level Path growing algorithm followed by dynamic programming
// timer[18][]	2 round suitor-based weighted matching + DP, using sorted edge lists
// timer[20][]	Sorting of edge list of each vertex
// timer[21][]	Greedy weighted algorithm, considering vertices in order, compare with [3]

// Parallel timings, all timings gives the best run
// ****************
// ** Unweighted algorithms **
// p_timer[0][i][conf]	Unweighted greedy vertex based matching algorithm
// p_timer[1][i][conf]	Unweighted edge based verification code, does not return matching
// p_timer[2][i][conf]	Unweighted edge based locking code
// p_timer[3][i][conf]	Not in use
// ** Weighted algorithms **
// p_timer[4][i][conf]	Parallel suitor algorithm
// p_timer[5][i][conf]	2 stage weighted suitor algorithm with dynamic programming
// p_timer[6][i][conf]	Parallel suitor algorithm with sorted adjacency lists
// p_timer[7][i][conf]	2 stage weighted suitor algorithm with sorted edge lists followed by dynamic programming
//
// p_cost match up with p_timer for 4,5,6, and 7 gives the cost of the matchings


  store_value(wf,"TimeSuitorSequential",timer[5],n_graphs);    // Print timing results for sequential suitor algorithm
  store_value(wf,"CostSuitorSequential",cost[5],n_graphs);     // Print cost results for sequential suitor algorithm

  store_value(wf,"TimeSortedSuitorSequential",timer[13],n_graphs); // Print timing results for sequential suitor algorithm with sorted edge lists
  store_value(wf,"CostSortedSuitorSequential",cost[13],n_graphs);  // Print cost results for sequential suitor algorithm with sorted edge lists

  store_value(wf,"Time2RSuitorSequential",timer[15],n_graphs);  // Print timing results for 2 round sequential suitor algorithm
  store_value(wf,"Cost2RSuitorSequential",cost[15],n_graphs);   // Print cost results for 2 round sequential suitor algorithm

  store_value(wf,"Time2RGreedySequential",timer[16],n_graphs);  // Print timing results for 2 round sequential suitor algorithm
  store_value(wf,"Cost2RGreedySequential",cost[16],n_graphs);   // Print cost results for 2 round sequential suitor algorithm

  store_value(wf,"TimeGPASequential",timer[11],n_graphs);      // Print timing results for sequential GPA algorithm
  store_value(wf,"CostGPASequential",cost[11],n_graphs);       // Print cost results for sequential GPA algorithm

  store_value(wf,"TimePathGrowSequential",timer[12],n_graphs); // Print timing results for sequential Path growing algorithm
  store_value(wf,"CostPathGrowSequential",cost[12],n_graphs);  // Print cost results for sequential Path growing algorithm

  store_value(wf,"TimeGreedySequential",timer[21],n_graphs); // Print timing results for sequential greedy algorithm
  store_value(wf,"CostGreedySequential",cost[21],n_graphs);  // Print cost results for sequential greedy algorithm

  store_value(wf,"TimeEdgeSortSequential",timer[10],n_graphs); // Print timing results for sorting all edges
  store_value(wf,"TimeEdgeListsSortSequential",timer[20],n_graphs); // Print timing results for sorting each edge list
  
/*
  fprintf(wf,"s_vx = [");       // Print results for sequential vertex based matching to file
  for(i=0;i<n_graphs;i++)  {
    fprintf(wf,"%lf ",timer[0][i]);
  }
  fprintf(wf,"];\n");

  fprintf(wf,"s_ed = [");       // Print results for sequential edge based matching to file
  for(i=0;i<n_graphs;i++)  {
    fprintf(wf,"%lf ",timer[1][i]);
  }
  fprintf(wf,"];\n");

  fprintf(wf,"p_vx = [");       // Print results for parallel vertex based matching to file
  for(i=0;i<n_graphs;i++)  {
    for(j=0;j<n_conf;j++)  {
      fprintf(wf,"%lf ",t_ver[i][j]);
    }
    if (i != n_graphs-1)
      fprintf(wf,";\n");
  }
  fprintf(wf,"];\n");

  fprintf(wf,"p_ex = [");       // Print results for parallel edge based matching to file
  for(i=0;i<n_graphs;i++)  {
    for(j=0;j<n_conf;j++)  {
      fprintf(wf,"%lf ",t_edge[i][j]);
    }
    if (i != n_graphs-1)
      fprintf(wf,";\n");
  }
  fprintf(wf,"];\n");

  fprintf(wf,"p_lc = [");       // Print results for parallel lock based matching to file
  for(i=0;i<n_graphs;i++)  {
    for(j=0;j<n_conf;j++)  {
      fprintf(wf,"%lf ",t_lock[i][j]);
    }
    if (i != n_graphs-1)
      fprintf(wf,";\n");
  }
  fprintf(wf,"];\n");
*/
  fclose(wf);
  printf("Normal exit\n");
}

int read_graph(int *n,int *m,char *f_name,int **ver,int **edges,edge **e,
               double **weight,wed **we,neig **swe) {

// Note reuse of p as count and left as start

  FILE *fp;
  long int i;
  int x,y,z;
  double v;
  MM_typecode matcode;
  int cols;

  int *count;
  int *start;

  // printf("Opening file %s \n",f_name);

// First make sure graph is of the right type

  fp = fopen(f_name,"r");
  if (fp == NULL) {
    printf("Could not open file %s \n",f_name);
    return(false);
  }

  if (mm_read_banner(fp, &matcode) != 0) {
    printf("Could not process Matrix Market banner.\n");
    return(false);
  }

  if ((mm_read_mtx_crd_size(fp, &cols, n, m)) !=0) {
    printf("Could not read size of graph.\n");
    return(false);
  }

  if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) || !mm_is_sparse(matcode) || !mm_is_symmetric(matcode)) {
    printf("The program can only read files that are sparse symmmetric matrices in coordinate format! \n");
    return(false);
  }

// fscanf(fp,"%d %d",&n,&m);
  printf("Reading graph named %s containing %d vertices and %d edges\n",f_name,*n,*m);

// Pointers into edge list for each vertex is stored in ver

  *ver = (int *) malloc(sizeof(int)*(2+(*n)));
  if (*ver == NULL) {
    printf("Unable to allocate space for vertex-array in read_graph \n");
    return(false);
  }

// 'edges' contains the compressed edge lists of each vertex, each edge is stored twice

  *edges = (int *) malloc(2*sizeof(int)*(*m));
  if (*edges == NULL) {
    printf("Unable to allocate space for edges-array in read_graph \n");
    return(false);
  }

// 'weight' contains the weight of each edge, stored in the same way as the edge-lists

  *weight= (double *) malloc(2*sizeof(double)*(*m));
  if (*weight== NULL) {
    printf("Unable to allocate space for weight-array in read_graph \n");
    return(false);
  }

// The raw edge list is stored in e, i.e. in the same order as it was read in

  *e = (edge *) malloc(sizeof(edge)*(*m));
  if (*e == NULL) {
    printf("Unable to allocate space for e-array in read_graph \n");
    return(false);
  }

// The edge list with the weight

  *we = (wed *) malloc((*m) * sizeof(wed));
  if (*we == NULL) {
    printf("Unable to allocate space for we-array in read_graph \n");
    return(false);
  }

// The compressed edge list with the weight, needed for sorting using qsort()

  *swe = (neig *) malloc(((*m)*2) * sizeof(neig));
  if (*swe == NULL) {
    printf("Unable to allocate space for swe-array in read_graph \n");
    return(false);
  }

// Pointers used for putting the edges in their correct places

  start = (int *) malloc(sizeof(int)*((*n)+1));
  if (start == NULL) {
    printf("Unable to allocate space for start-array in read_graph \n");
    return(false);
  }

// Used for calculating where the edge list of each vertex should be stored in the compressed edge list

  count = (int *) malloc(sizeof(int)*((*n)+1));
  if (count == NULL) {
    printf("Unable to allocate space for count-array in read_graph \n");
    return(false);
  }


// Start counting of degrees by setting counters to zero
  for(i=1;i<=*n;i++) {
    count[i] = 0;
  }

  int num_edges = 0;

  printf("Number of possible edges is %d \n",*m);
  
// Read inn the edges
  if (mm_is_real(matcode)) {
    printf("Starting to read %d edges \n",*m);
    for(i=0;i<*m;i++) {
      fscanf(fp,"%d %d %lf",&x,&y,&v); // Use this line if there is a weight 
      if (x != y) { // Avoid self-edges
        (*e)[num_edges].x = x; // Store edges
        (*e)[num_edges].y = y;

        (*we)[num_edges].x = x;
        (*we)[num_edges].y = y;
        int intv = (int) fabs(v);
        (*we)[num_edges].w = (double) intv; 
//        if ((fabs(v) < 0.0) || (fabs(v) > 5793)) {
        if (fabs(v) < 0.0)  {
          printf("error reading data,got %f \n",fabs(v));
          return;
        }
      //   (*we)[num_edges].w = 1.0;
      //  (*we)[num_edges].w = (double)rand();

        count[x]++; // Get vertex degrees
        count[y]++;

        num_edges++;
        if (num_edges > *m) {
          printf("Have set num_edges to %d while i=%d \n",num_edges,i);
          return;
        }
      }
    }
  }
  else {          // Symbolic matrix
    printf("Trouble ahead, the code now assumes weighted graphs \n");
    for(i=0;i<*m;i++) {
      fscanf(fp,"%d %d",&x,&y);
      if (x != y) { // Avoid self-edges
        (*e)[num_edges].x = x; // Store edges
        (*e)[num_edges].y = y;

        count[x]++; // Get vertex degrees
        count[y]++;

        num_edges++;
      }
    }
  }
  printf("original edges %d, used edges %d \n",*m,num_edges);

  *m = num_edges;  // Make sure m is the correct number of edges

// Find starting positions in edge list for each vertex
  start[1] = 0;
  (*ver)[1] = 0;
  for(i=2;i<=*n+1;i++) {
    start[i] = start[i-1]+count[i-1];
    (*ver)[i] = start[i];
  }

// Place edges in edge lists, once for each endpoint

  for(i=0;i<*m;i++) {
    x = (*we)[i].x;
    y = (*we)[i].y;
    v = (*we)[i].w;
    if ((x == 0) || (y == 0)) {
      printf("edge %d: %d and %d are neighbors, numbering starts from 1! \n",i,x,y);
      return(false);
    }

    (*edges)[start[x]] = y;
    (*weight)[start[x]] = v;
    start[x]++;

    (*edges)[start[y]] = x;
    (*weight)[start[y]] = v;
    start[y]++;

  }

  fclose(fp);
  free(count);
  free(start);
  return(true);
}

// Get input. This is given as a file name to the executable. 
// This file should contain the following information (per line):
// Number of input graphs
// Number of times each configuration is run, the program reports the best running time out of these
// Number of configurations, followed by the number of threads in each configuration
// One line giving the complete name of each graph file that is to be run

int get_input(int argc, char *argv[],int *n_graphs,int *n_runs,int *n_conf,int **conf,char ***name) {

  FILE *rf;	// File pointer
  int i;

  if (argc != 2) {
    printf("Give data file name as first input parameter!\n");
    return(false);
  }

  printf("Reading setup from %s \n",argv[1]);

// Opening file containing file names for graphs

  rf = fopen(argv[1],"r");
  if (rf == NULL) {
    printf("Cannot open data file: %s \n",argv[1]);
    return(false);
  }


// Get number of graphs to read
  fscanf(rf,"%d",n_graphs);

  *name =  malloc(sizeof(char *)*(*n_graphs));  // Allocate one pointer for each graph 

  if (*name == NULL) {
    printf("Unable to allocate space for names of %d graphs\n",max_graphs);
    return(0);
  }

// Get number of runs per configuration
  fscanf(rf,"%d",n_runs);

// Get number of thread configurations
  fscanf(rf,"%d",n_conf);
  *conf = (int *) malloc(sizeof(int)*(*n_conf));  // Allocate space for configurations


  if (*conf == NULL) {
    printf("Unable to allocate memory for %d different thread configurations \n",n_conf);
    return(0);
  }

// Get the different configurations
  for(i=0;i< *n_conf;i++) {
    fscanf(rf,"%d",&(*conf)[i]);
  }

// Get the different file names
  for(i=0;i< *n_graphs;i++) {
    (*name)[i] = (char *) malloc(sizeof(char)*100);
    if ((*name)[i] == NULL) {
      printf("Unable to allocate memory for graph name %d \n",i);
      return(0);
    }
// Read name of graph
    fscanf(rf,"%s",(*name)[i]);
//    printf("Read name %s \n",(*name)[i]);
  }

  fclose(rf);
  return(true);
}


int allocate_memory(int size,int **n,int **m,double ***timer,double ***cost,double ****p_timer,double ****p_cost,int n_conf) {

  *n = (int *) malloc(sizeof(int)*size);  // List holding the number of vertices in each graph

  if (*n == NULL) {
    printf("Unable to allocate memory for n[] in allocate_memory() \n");
    return(false);
  }

  *m = (int *) malloc(sizeof(int)*size);  // List holding the number of edges in each graph
  if (*m == NULL) {
    printf("Unable to allocate memory for m[] in allocate_memory() \n");
    return(false);
  }

  *timer = malloc(sizeof(double *)*(max_experiment));  // Allocate one pointer for each sequential experiment
  if (*timer == NULL) {
    printf("Unable to allocate memory for timer in allocate_memory() \n");
    return(false);
  }

  *cost = malloc(sizeof(double *)*(max_experiment));  // Allocate one pointer for each sequential experiment
  if (*timer == NULL) {
    printf("Unable to allocate memory for cost in allocate_memory() \n");
    return(false);
  }

  // For each experiment, allocate one number for each graph
  int i;
  for(i=0;i<max_experiment;i++) {
    (*cost)[i] = (double *) malloc(sizeof(double)*size);
    if ((*cost)[i] == NULL) {
      printf("Unable to allocate memory for cost %d in allocate_memory() \n",i);
      return(false);
    }
    (*timer)[i] = (double *) malloc(sizeof(double)*size);
    if ((*timer)[i] == NULL) {
      printf("Unable to allocate memory for timer %d in allocate_memory() \n",i);
      return(false);
    }
    int j;
    for(j=0;j<size;j++)          // Set each timer to -1
      (*timer)[i][j] = -1.0;
  }

  *p_timer = malloc(sizeof(double **)*(max_experiment));  // Allocate one pointer for each parallel experiment
  *p_cost  = malloc(sizeof(double **)*(max_experiment));  // Allocate one pointer for each parallel experiment

  if (*p_timer == NULL) {
    printf("Unable to allocate memory for p_timer in allocate_memory() \n");
    return(false);
  }
  if (*p_cost == NULL) {
    printf("Unable to allocate memory for p_cost in allocate_memory() \n");
    return(false);
  }
  // For each experiment, allocate "size" pointers for each graph
  for(i=0;i<max_experiment;i++) {
    (*p_timer)[i] = (double **) malloc(sizeof(double *)*size);
    (*p_cost)[i]  = (double **) malloc(sizeof(double *)*size);
    if ((*p_timer)[i] == NULL) {
      printf("Unable to allocate memory for p_timer %d in allocate_memory() \n",i);
      return(false);
    }
    if ((*p_cost)[i] == NULL) {
      printf("Unable to allocate memory for p_cost %d in allocate_memory() \n",i);
      return(false);
    }
    int j;
    for(j=0;j<size;j++) {
      (*p_timer)[i][j] = (double *)  malloc(sizeof(double)*n_conf);
      (*p_cost)[i][j] = (double *)  malloc(sizeof(double)*n_conf);
      if ((*p_timer)[i][j] == NULL) {
        printf("Unable to allocate memory for p_timer %d in allocate_memory(),%d \n",i,j);
        return(false);
      }
      if ((*p_cost)[i][j] == NULL) {
        printf("Unable to allocate memory for p_cost %d,%d in allocate memory \n",i,j);
        return(false);
      }
      int k;
      for(k=0;k<n_conf;k++) {
        (*p_timer)[i][j][k] = -1.0;
        (*p_cost)[i][j][k] = -1.0;
      }
    }
  }
  return(true);
}


int prepare_output(FILE **wf,int n_conf,int conf[]) {

// Open file for writing of results
  *wf = fopen("results.m","w");
  if (*wf == NULL) {
    printf("Unable to open results.m for writing \n");
    return(false);
  }

// print the thread configurations to file
  fprintf(*wf,"x = [");
  int i;
  for(i=0;i<n_conf;i++) 
    fprintf(*wf,"%d ",conf[i]);
  fprintf(*wf,"];\n");

  return(true);
}

// Verify that p[] defines a legal matching

int verify_matching(int n,int *ver,int *edges,int *p) {

  int i;

  for(i=1;i<=n;i++) {
    if ((p[i] < 0) || (p[i] > n)) {
      printf("p[%d] = %d, while n = %d \n",i,p[i],n);
      return(false);
    }
    if ((p[i] != 0) && (p[p[i]] != i)) {
      printf("p[%d] = %d, while p[%d] = %d \n",i,p[i],p[i],p[p[i]]); 
      return(false);
    }
  }
  return(true);
}

// Compute a random ordering of the vertices

int random_order(int n,int *order) {

  int l;

// Compute a random ordering of the vertices
  for(l=1;l<=n;l++) {
    order[l] = l;
  }

  srandom(1);
  for(l=1;l<n;l++) {
    long int x = random() % (long int) (n-l+1);
    x++;
    int tmp = order[n-l+1];
    order[n-l+1] = order[x];
    order[x] = tmp;
  }
  return(true);
}

// Allocating memory specifically for this graph

int allocate_graph_memory(int n,int **next,int **used,int **p,int **p1,int **p2) {		

  *p = (int *) malloc(sizeof(int)*(n+2));	// Storage for final matching
  if (*p == NULL) {
    printf("Unable to allocate space for p[] in allocate_graph_memory() \n");
    return(false);
  }

  *p1 = (int *) malloc(sizeof(int)*(n+2));	// Storage for first matching in two level algorithms
  if (*p1 == NULL) {
    printf("Unable to allocate space for p1[] in allocate_graph_memory() \n");
    return(false);
  }

  *p2 = (int *) malloc(sizeof(int)*(n+2));	// Storage for second matching in two level algorithms
  if (*p2 == NULL) {
    printf("Unable to allocate space for p2[] in allocate_graph_memory() \n");
    return(false);
  }

  *next = (int *) malloc(sizeof(int)*n*2);	// Pointers used in suitor based algorithm with sorted edge lists
  if (*next == NULL) {
    printf("Unable to allocate space for next[] in allocate_graph_memory() \n");
    return(false);
  }

  *used = (int *) malloc(sizeof(int)*(n+2));
  if (*used == NULL) {
    printf("Unable to allocate space for used[] in allocate_graph_memory() \n");
    return(false);
  }
  return(true);
}


// Calculate the cost of the matching

double cost_matching(int n,int *ver,int *edges,double *weight,int *match) {

  double glob_sum = 0.0;
  int i,k;

  for(i=1;i<=n;i++) {
    if ((match[i] != 0) && (i < match[i])) { // Only use vertices that are matched and that have lower index then their partner
      for(k=ver[i];k<ver[i+1];k++) {         // Loop through neighbors of vertex i
        if (edges[k] == match[i]) {
          glob_sum += weight[k];
          if (match[edges[k]] != i)
            printf("Error in cost_matching: %d is matched with %d but %d is matched with %d \n",i,match[i],edges[k],match[edges[k]]);
        } // if
      } // for k
    } // if
  } // for i
  return(glob_sum);
}

void store_value(FILE *wf,char *s,double *values,int n) {
  int i;
  fprintf(wf,"%s = [",s); 
  for(i=0;i<n;i++)  {
    fprintf(wf,"%lf ",values[i]);
  }
  fprintf(wf,"];\n");
}
