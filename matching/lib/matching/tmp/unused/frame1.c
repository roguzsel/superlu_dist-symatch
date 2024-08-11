#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "frame1.h"
#include "seq1.c"
#include "eseq1.c"
#include "ver1.c"
#include "edge1.c"
#include "lock1.c"
#include "weight1.c"
#include "mmio.c"

  int p[max_node];            // If p(i) = j then i thinks it is matched with j
  int left[max_node];                  // Storage for vertices that are to be reexamined

int main(int argc, char *argv[]) {


  edge *e;           // List of all edges
  int *ver;                   // Pointers to edge lists
  int *edges;                 // Edge lists

  int m[max_graphs];                      // Number of edges
  int n[max_graphs];                      // Number of vertices

  int read_graph();
  void discontinue();
  int i,j;

  FILE *rf;
  FILE *wf;
  char name[max_graphs][100];   // Name of current graph
  char *tc;         // Name of current graph without .mtx
  int n_graphs;     // Number of graphs to run
  int n_runs;       // Number of runs for each configuration
  int n_conf;       // Number of thread configurations
  int conf[100];    // Number of threads in each configuration

  float t_bs[max_graphs];          // timings for sequential vertex oriented algorithm
  float t_bse[max_graphs];         // timings for sequential edge based algorithm
  float t_edge[max_graphs][100];   // timings for parallel edge based algorithm 
  float t_ver[max_graphs][100];    // timings for parallel vertex based algorithm
  float t_lock[max_graphs][100];   // timings for parallel edge based locking algorithm
  float t_weight[max_graphs][100]; // timings for parallel weighted algorithm

  double mt1,mt2;

  printf("Starting algorithm \n");

  if (argc != 2) {
    discontinue("Give data file name as parameter\n");
  }

  printf("starting reading from %s \n",argv[1]);

// Opening file containing file names for graphs
  rf = fopen(argv[1],"r");
  if (rf == NULL)
    discontinue("Cannot open data file \n");

// Get number of graphs to read
  fscanf(rf,"%d",&n_graphs);
  if (n_graphs > max_graphs) {
    printf("Can only handle %d graphs\n",max_graphs);
    return (0);
  }

// Get number of runs per configuration
  fscanf(rf,"%d",&n_runs);
// Get number of thread configurations
  fscanf(rf,"%d",&n_conf);
// Get the different configurations
  for(i=0;i<n_conf;i++) {
    fscanf(rf,"%d",&(conf[i]));
  }

// Open file for writing of results
  wf = fopen("results.m","w");
  if (wf == NULL) {
    printf("Unable to open results for writing \n");
  }

// print the thread configurations to file
  fprintf(wf,"x = [");
  int jj;
  for(jj=0;jj<n_conf;jj++) 
    fprintf(wf,"%d ",conf[jj]);
  fprintf(wf,"];\n");

// Run the main loop 
  for(i=0;i<n_graphs;i++) {
    t_bs[i] = -1;
    t_bse[i] = -1;

// Read name of graph
    fscanf(rf,"%s",name[i]);

    printf("Starting to read graph number %d: %s \n",i,name[i]);

// Reading graph
    if (!read_graph(&(n[i]),&(m[i]),name[i],&ver,&edges,&e)) {
      printf("Problem reading graph %s \n",name[i]);
    }
/*
  left = (int *) malloc(sizeof(int)*n[0]);
  if (left == NULL) {
    printf("Unable to allocate space for left-array \n");
    return(FALSE);
  }
*/
    tc = strtok(name[i],".");

// Run sequential code 
    int l;
    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      seq(n[i],ver,edges,p);
      mt2 = omp_get_wtime(); 
      if ((t_bs[i] < 0) || (mt2-mt1 < t_bs[i]))
        t_bs[i] = mt2-mt1;
    }
    printf("Sequential vertex based greedy matching took %f seconds\n",t_bs[i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      eseq(n[i],m[i],e,p);
      mt2 = omp_get_wtime(); 
      if ((t_bse[i] < 0) || (mt2-mt1 < t_bse[i]))
        t_bse[i] = mt2-mt1;
    }

    printf("Sequential edge based greedy matching took %f seconds\n",t_bse[i]);

// Allocate space for locks

    omp_lock_t *nlocks;
    nlocks = (omp_lock_t *) malloc(n[i]*sizeof(omp_lock_t));

// Run parallel algorithm
// One iteration for each thread configuration
    for(j=0;j<n_conf;j++) {
      printf("Starting configuration %d \n",j);
      omp_set_num_threads(conf[j]);  // Set number of threads in this configuration

// Initialize timing variables

      t_edge[i][j] = -1;
      t_ver[i][j] = -1;
      t_lock[i][j] = -1;
      t_weight[i][j] = -1;

    

#pragma omp parallel 
      {
        int k;
        int threads = omp_get_num_threads();
        int my_id = omp_get_thread_num();
      
      // int *tp = (int *) malloc(n[i]*sizeof(int));


        if (threads != conf[j])
          printf("***** Did not get the correct numer of threads! Wanted %d, got %d *****\n",conf[i],threads);

// Now we are ready to run the different algorithms
// Run as many times as required, only store best run

        for(k=0;k<n_runs;k++) {
// Running vertex based verification code

#pragma omp barrier
#pragma omp master
          {  mt1 = omp_get_wtime(); }
          verify(n[i],ver,edges,p,left);
#pragma omp master
          {  mt2 = omp_get_wtime();
             if ((t_ver[i][j] < 0) || (mt2-mt1 < t_ver[i][j]))
               t_ver[i][j] = mt2-mt1;
          }
#pragma omp barrier

// Running edge based verification code

#pragma omp master
          {  mt1 = omp_get_wtime(); }
             edgec(n[i],m[i],e,p,left);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t_edge[i][j] < 0) || (mt2-mt1 < t_edge[i][j]))
               t_edge[i][j] = mt2-mt1;
          }

// Running edge based locking code

#pragma omp master
          {  mt1 = omp_get_wtime(); }
             locking1(n[i],m[i],e,p,nlocks);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t_lock[i][j] < 0) || (mt2-mt1 < t_lock[i][j]))
               t_lock[i][j] = mt2-mt1;
          }

// Running weighted algorithm

#pragma omp master
          {  mt1 = omp_get_wtime(); }
           //  weighted(n[i],ver,edges,p,left);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t_weight[i][j] < 0) || (mt2-mt1 < t_weight[i][j]))
               t_weight[i][j] = mt2-mt1;
          }

        } // End of k iterations over same configuration 
       
#pragma omp master
        { printf("Vertex code using %d threads took %f seconds, speedup = %4.2f \n",conf[j],t_ver[i][j],t_bs[i]/t_ver[i][j]);
          printf("Edge code using %d threads took %f seconds, speedup = %4.2f \n",conf[j],t_edge[i][j],t_bse[i]/t_edge[i][j]);
          printf("Locking code using %d threads took %f seconds, speedup = %4.2f \n",conf[j],t_lock[i][j],t_bse[i]/t_lock[i][j]);
          printf("One sided code using %d threads took %f seconds, speedup = %4.2f \n",conf[j],t_weight[i][j],t_bse[i]/t_weight[i][j]);
        }

      }  // End of parallel section


    }  // End of different thread configurations

    printf("Freeing up space \n");
    free(nlocks);   // Free up space used for locks
    free(ver);
    free(edges);
    free(e);
    printf("Done freeing space\n");

  } // Loop over graphs

// Print results to file

  printf("Printing to file \n");
  fprintf(wf,"name = char(");       // Print the names of the graphs
  for(i=0;i<n_graphs;i++)  {
    fprintf(wf,"'%s' ",name[i]);
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

  fprintf(wf,"s_vx = [");       // Print results for sequential vertex based matching to file
  for(i=0;i<n_graphs;i++)  {
    fprintf(wf,"%f ",t_bs[i]);
  }
  fprintf(wf,"];\n");

  fprintf(wf,"s_ed = [");       // Print results for sequential edge based matching to file
  for(i=0;i<n_graphs;i++)  {
    fprintf(wf,"%f ",t_bse[i]);
  }
  fprintf(wf,"];\n");

  fprintf(wf,"p_vx = [");       // Print results for parallel vertex based matching to file
  for(i=0;i<n_graphs;i++)  {
    for(j=0;j<n_conf;j++)  {
      fprintf(wf,"%f ",t_ver[i][j]);
    }
    if (i != n_graphs-1)
      fprintf(wf,";\n");
  }
  fprintf(wf,"];\n");

  fprintf(wf,"p_ex = [");       // Print results for parallel edge based matching to file
  for(i=0;i<n_graphs;i++)  {
    for(j=0;j<n_conf;j++)  {
      fprintf(wf,"%f ",t_edge[i][j]);
    }
    if (i != n_graphs-1)
      fprintf(wf,";\n");
  }
  fprintf(wf,"];\n");

  fprintf(wf,"p_lc = [");       // Print results for parallel lock based matching to file
  for(i=0;i<n_graphs;i++)  {
    for(j=0;j<n_conf;j++)  {
      fprintf(wf,"%f ",t_lock[i][j]);
    }
    if (i != n_graphs-1)
      fprintf(wf,";\n");
  }
  fprintf(wf,"];\n");

  fclose(wf);
  fclose(rf);
  printf("Normal exit\n");
}

int read_graph(int *n,int *m,char *f_name,int **ver,int **edges,edge **e) {

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

// 'edges' contains the compressed edge lists of each vertex

  *edges = (int *) malloc(sizeof(int)*(*m));
  if (*edges == NULL) {
    printf("Unable to allocate space for edges-array in read_graph \n");
    return(FALSE);
  }

// The raw edge list is stored in e, i.e. in the same order as it was read in

  *e = (edge *) malloc(sizeof(edge)*(*m));
  if (*e == NULL) {
    printf("Unable to allocate space for e-array in read_graph \n");
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

        if ((num_edges % 2) == 0)
          count[x]++; // Get vertex degrees
        else
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

        if ((num_edges % 2) == 0)
          count[x]++; // Get vertex degrees
        else
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
    if ((x == 0) || (y == 0)) {
      printf("%d and %d er naboer \n",x,y);
      return(FALSE);
    }
    if ((i % 2) == 0) {
      (*edges)[start[x]] = y;
      start[x]++;
    }
    else {
      (*edges)[start[y]] = x;
      start[y]++;
    }
  }

  fclose(fp);
  return(TRUE);
}

void discontinue(char* s)  {
  printf("%s",s);
  exit(1);
}
