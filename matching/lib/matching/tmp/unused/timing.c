#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "frame1.h"
#include "mmio.c"



#define big_number 10000000
#define small_number 100
#define rep 5

int main(int argc, char *argv[]) {
  FILE *rf;
  int i;
  int p=1;
  int result[small_number];
  float t[100];
  float mt1,mt2;
//  int value[big_number];
  int *value;
  int r,x;

  int read_graph();

  int n,m;
  char name[max_graphs][100];   // Name of current graph
  int n_graphs;
  int n_runs;
  int n_conf;
  int conf[100];    // Number of threads in each configuration
  int *ver;                   // Pointers to edge lists
  int *edges;                 // Edge lists
  float *weight;              // Corresponding edge weights
  edge *e;           // List of all edges
  float *w;          // Corresponding weight of each edge


  printf("Starting algorithm \n");

  if (argc != 2) {
    printf("Give data file name as parameter\n");
    return(0);
  }

  printf("starting reading from %s \n",argv[1]);

// Opening file containing file names for graphs
  rf = fopen(argv[1],"r");
  if (rf == NULL) {
    printf("Cannot open data file \n");
    return(0);
  }

// Get number of graphs to read
  fscanf(rf,"%d",&n_graphs);
  if (n_graphs > max_graphs) {
    printf("Can only handle %d graphs\n",max_graphs);
    return(0);
  }

// Get number of runs per configuration
  fscanf(rf,"%d",&n_runs);
// Get number of thread configurations
  fscanf(rf,"%d",&n_conf);
// Get the different configurations
  for(i=0;i<n_conf;i++) {
    fscanf(rf,"%d",&(conf[i]));
  }

// Read name of graph
    fscanf(rf,"%s",name[0]);

    printf("Starting to read graph number %d: %s \n",i,name[0]);


  if (!read_graph(&n,&m,name[0],&ver,&edges,&e,&w,&weight)) {
    printf("Problem reading graph %s \n",name[0]);
  }

  value = (int *) malloc(big_number*sizeof(int));

  for(i=0;i<n_conf;i++) {
     omp_set_num_threads(conf[i]);
     result[i] = 0;
     t[i] = -1;

     for(r=0;r<n_runs;r++) {
#pragma omp parallel 
      {
        int k,j;
        int threads = omp_get_num_threads();
        int my_id = omp_get_thread_num();
/*
        if (my_id == 0) 
          printf("Using %d threads\n",threads);
*/

#pragma omp barrier

#pragma omp master
        { mt1 = omp_get_wtime(); }

// Do work

        int sum;
#pragma omp for private(k,j,sum)
        for(k=0;k<n;k++) {
          sum = 0;
          for(j=ver[k];j<ver[k+1];j++) {
            sum = (edges[j]) % 13;
            value[k] += sum;
          }
        }

#pragma omp master
        {  mt2 = omp_get_wtime();
           if ((mt2-mt1 < t[i]) || (t[i] < 0))
             t[i] = mt2-mt1;
           for(k=0;k<small_number;k++) {
             result[i] += value[k];
           }
        }
        #pragma omp barrier


      }  // End of parallel section
    }  // End of number of repetitions
  }  // End of configurations

  for(i=0;i<n_conf;i++) {
    printf("%d: Time %f, Speedup %f \n",conf[i],t[i],t[0]/t[i]);
  }
  int sum=0;
  int k;
  for(k=0;k<small_number;k++) 
    sum += result[k];
  printf("%d \n",sum);
}

int read_graph(int *n,int *m,char *f_name,int **ver,int **edges,edge **e,float **w,float **weight) {

// Note reuse of p as count and left as start

  FILE *fp;
  long int i;
  int x,y,z;
  float v;
  MM_typecode matcode;
  int cols;
  int *count;
  int *start;

  printf("Opening file %s \n",f_name);

  fp = fopen(f_name,"r");
  if (fp == NULL) {
    printf("Could not open file %s \n",f_name);
    return(FALSE);
  }

  if (mm_read_banner(fp, &matcode) != 0) {
    printf("Could not process Matrix Market banner.\n");
    return(FALSE);
  }

  if ((mm_read_mtx_crd_size(fp, &cols, n, m)) !=0) {
    printf("Could not read size of graph.\n");
    return(FALSE);
  }

  if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) || !mm_is_sparse(matcode) || !mm_is_symmetric(matcode)) {
    printf("The program can only read files that are sparse symmmetric matrices in coordinate format! \n");
    return(FALSE);
  }

// fscanf(fp,"%d %d",&n,&m);
  printf("Reading graph named %s containing %d vertices and %d edges\n",f_name,*n,*m);

// Pointers into edge list for each vertex is stored in ver

  *ver = (int *) malloc(sizeof(int)*(2+(*n)));
  if (*ver == NULL) {
    printf("Unable to allocate space for vertex-array in read_graph \n");
    return(FALSE);
  }

// 'edges' contains the compressed edge lists of each vertex, each edge is stored twice

  *edges = (int *) malloc(2*sizeof(int)*(*m));
  if (*edges == NULL) {
    printf("Unable to allocate space for edges-array in read_graph \n");
    return(FALSE);
  }

// 'weight' contains the weight of each edge, stored in the same way as the edge-lists

  *weight= (float *) malloc(2*sizeof(float)*(*m));
  if (*weight== NULL) {
    printf("Unable to allocate space for weight-array in read_graph \n");
    return(FALSE);
  }

// The raw edge list is stored in e, i.e. in the same order as it was read in

  *e = (edge *) malloc(sizeof(edge)*(*m));
  if (*e == NULL) {
    printf("Unable to allocate space for e-array in read_graph \n");
    return(FALSE);
  }

// The raw weight list is stored in w, i.e. in the same order as it was read in

  *w = (float *) malloc(sizeof(float)*(*m));
  if (*w == NULL) {
    printf("Unable to allocate space for w-array in read_graph \n");
    return(FALSE);
  }

// Pointers used for putting the edges in their correct places

  start = (int *) malloc(sizeof(int)*(*n));
  if (start == NULL) {
    printf("Unable to allocate space for start-array in read_graph \n");
    return(FALSE);
  }

// Used for calculating where the edge list of each vertex should be stored in the compressed edge list

  count = (int *) malloc(sizeof(int)*(*n));
  if (count == NULL) {
    printf("Unable to allocate space for count-array in read_graph \n");
    return(FALSE);
  }

// Don't remember why this is here...

  if ((*n>max_node) || (*m > max_edge)) {
    printf("Maximum number of nodes and edges is %d and %d \n",max_node,max_edge);
    return(FALSE);
  }

// Start counting of degrees by setting counters to zero
  for(i=1;i<=*n;i++) {
    count[i] = 0;
  }

  int num_edges = 0;

// Read inn the edges
  if (mm_is_real(matcode))
    for(i=0;i<*m;i++) {
      fscanf(fp,"%d %d %f",&x,&y,&v); // Use this line if there is a weight 
      if (x != y) { // Avoid self-edges
        (*e)[num_edges].x = x; // Store edges
        (*e)[num_edges].y = y;
        (*w)[num_edges] = abs(v);

//        if ((num_edges % 2) == 0)
        count[x]++; // Get vertex degrees
//        else
        count[y]++;

        num_edges++;
      }
    }
  else           // Symbolic matrix
    for(i=0;i<*m;i++) {
      fscanf(fp,"%d %d",&x,&y);
      if (x != y) { // Avoid self-edges
        (*e)[num_edges].x = x; // Store edges
        (*e)[num_edges].y = y;

//        if ((num_edges % 2) == 0)
        count[x]++; // Get vertex degrees
//        else
        count[y]++;

        num_edges++;
      }
    }

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
    x = (*e)[i].x;
    y = (*e)[i].y;
    v = (*w)[i];
    if ((x == 0) || (y == 0)) {
      printf("%d and %d er naboer \n",x,y);
      return(FALSE);
    }
//    if ((i % 2) == 0) {
    (*edges)[start[x]] = y;
    (*weight)[start[x]] = v;
    start[x]++;
//    }
//    else {
    (*edges)[start[y]] = x;
    (*weight)[start[y]] = v;
    start[y]++;
//    }
  }

  fclose(fp);
  return(TRUE);
}

