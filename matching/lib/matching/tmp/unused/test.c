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
#include "localmax.c" 		// Round based local max algorithm
#include "plocalmax.c" 		// Parallel round based local max algorithm, using a global edge list
#include "plocalmax1.c" 	// Parallel round based local max algorithm, using local edge lists

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
    if  (v1->w < v2->w)
      return true;
    if  (v1->w > v2->w)
      return false;
    return (int) (v1->x < v2->x);
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

//  omp_lock_t *nlocks;
  int tlist[max_threads];  // List for storing data for each thread (in plocalmax)

  int read_graph(),get_input(),allocate_memory();
  int read_bin_graph();
  double cost_matching();
  void store_value();
  void store_p_value();

  int i,j;          // Loop indices
  double mt1,mt2;   // Used for timing of each algorithm
  int xx;

  ws1 = (double *) malloc(sizeof(double)*(max_n+1));
  ws2 = (double *) malloc(sizeof(double)*(max_n+1));

  for(i=0;i<n_graphs;i++) {
    int l;
    wed *edgel;
    int xxx = ver[n[i]+1];
    omp_lock_t *nlocks;

  } // Loop over graphs
  
}




// Allocating memory specifically for this graph

int allocate_graph_memory(neig **swe,double **weight,edge **e,wed **we,int **ver,int **edges,int n,int **next,int **used,int **p,int **p1,int **p2) {		

// The compressed edge list with the weight, needed for sorting using qsort()

  long temp_long = 2;
  temp_long = temp_long * max_m * sizeof(neig);
/*
  *swe = (neig *) malloc(temp_long);
  if (*swe == NULL) {
    printf("Unable to allocate space for swe-array in allocate_graph_memory() \n");
    return(false);
  }
*/

// The edge list with the weight
//  printf("MAX_M = %d \n",max_m);

  *we = (wed *) malloc(max_m * sizeof(wed));
  if (*we == NULL) {
    printf("Unable to allocate space for we-array in allocate_graph_memory() \n");
    return(false);
  }


// 'weight' contains the weight of each edge, stored in the same way as the edge-lists
  *weight= (double *) malloc(2*sizeof(double)*(max_m));
  if (*weight== NULL) {
    printf("Unable to allocate space for weight-array in allocate_graph_memory() \n");
    return(false);
  }
// The raw edge list is stored in e, i.e. in the same order as it was read in


  *e = (edge *) malloc(sizeof(edge)*(max_m));
  if (*e == NULL) {
    printf("Unable to allocate space for e-array in allocate_graph_memory() \n");
    return(false);
  }


  *ver = (int *) malloc(sizeof(int)*(2+n));
  if (*ver == NULL) {
    printf("Unable to allocate space for vertex-array in allocate_graph_memory() \n");
    return(false);
  }

  *edges = (int *) malloc(2*sizeof(int)*(max_m));
  if (*edges == NULL) {
    printf("Unable to allocate space for edges-array in allocate_graph_memory \n");
    return(false);
  }

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
void store_p_value(FILE *wf,char *s,double **values,int n,int n_conf) {
  int i,j;
  fprintf(wf,"%s = [",s); 
  for(i=0;i<n;i++) {
    for(j=0;j<n_conf;j++)  {
      fprintf(wf,"%lf ",values[i][j]);
    }
  }
  fprintf(wf,"];\n");
}

void store_value(FILE *wf,char *s,double *values,int n) {
  int i;
  fprintf(wf,"%s = [",s); 
  for(i=0;i<n;i++)  {
    fprintf(wf,"%lf ",values[i]);
  }
  fprintf(wf,"];\n");
}
