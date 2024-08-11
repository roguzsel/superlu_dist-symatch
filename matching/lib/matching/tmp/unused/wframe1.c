#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "frame1.h"
#include "seq1.c"
#include "wseq1.c"
#include "wseq2.c"
#include "wseq3.c"
#include "eseq1.c"
#include "ver1.c"
#include "edge1.c"
#include "edge2.c"
#include "lock1.c"
#include "weight1.c"
#include "weight2.c"
#include "weight3.c"
#include "sweight.c"
#include "sweight1.c"
#include "sweight3.c"
#include "sweight4.c"
#include "sweight5.c"
#include "sweight6.c"
#include "wsort.c"
#include "mmio.c"
#include "gpa.c"
#include "path.c"

#define true  1
#define false 0

  int p[max_node];            // If p(i) = j then i thinks it is matched with j
  int p2[max_node];            // If p(i) = j then i thinks it is matched with j, round 2

  static int compedge(const void *m1, const void *m2) {
    wed *v1 = (wed *) m1; 
    wed *v2 = (wed *) m2; 
    return (int) (v1->w < v2->w);
  }


int main(int argc, char *argv[]) {

  edge *e;           // List of all edges
  float *w;          // Corresponding weight of each edge
  int *ver;                   // Pointers to edge lists
  int *edges;                 // Edge lists
  float *weight;              // Corresponding edge weights
  int *left;                  // Storage for vertices that are to be reexamined
  int *next;                  // Pointer for weighted matching
  float *ws;                    // Best suitor in weighted case
  float *ws2;                   // Best suitor in weighted case, round 2

  int m[max_graphs];                      // Number of edges
  int n[max_graphs];                      // Number of vertices

  int read_graph();
  int i,j;

  FILE *rf;
  FILE *wf;
  char name[max_graphs][100];   // Name of current graph
  char *tc;         // Name of current graph without .mtx
  int n_graphs;     // Number of graphs to run
  int n_runs;       // Number of runs for each configuration
  int n_conf;       // Number of thread configurations
  int conf[100];    // Number of threads in each configuration
  int *dist_e;                 // Distributed edge lists
  int *order;       // Stores a random permutation of the vertices
  int *used;       // Keeps track of which vertices has been used in the dynamic programming

  double t_bs[max_graphs];          // timings for sequential vertex oriented algorithm
  double t_bse[max_graphs];         // timings for sequential edge based algorithm
  double t_bsw[max_graphs];         // timings for sequential weighted algorithm
  double t_bgw[max_graphs];         // timings for greedy sequential weighted algorithm
  double t_rgw[max_graphs];         // timings for greedy sequential weighted algorithm
  double t_dgw[max_graphs];         // timings for greedy sequential weighted algorithm
  double t_5gw[max_graphs];         // timings for greedy sequential weighted algorithm
  double t_6gw[max_graphs];         // timings for greedy sequential weighted algorithm
  double t_gpa[max_graphs];         // timings for greedy sequential weighted algorithm
  double t_pat[max_graphs];         // timings for greedy sequential weighted algorithm

  double t_edge[max_graphs][100];   // timings for parallel edge based algorithm 
  double t1_edge[max_graphs][100];   // timings for parallel edge based algorithm 
  double t_ver[max_graphs][100];    // timings for parallel vertex based algorithm
  double t_lock[max_graphs][100];   // timings for parallel edge based locking algorithm
  double t_weight[max_graphs][100]; // timings for parallel weighted algorithm
  double t1_weight[max_graphs][100]; // timings for parallel weighted algorithm
  double t2_weight[max_graphs][100]; // timings for parallel weighted algorithm

  double mt1,mt2;

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
    t_bsw[i] = -1;
    t_bgw[i] = -1;
    t_rgw[i] = -1;
    t_dgw[i] = -1;
    t_5gw[i] = -1;
    t_6gw[i] = -1;
    t_gpa[i] = -1;
    t_pat[i] = -1;

// Read name of graph
    fscanf(rf,"%s",name[i]);

    printf("Starting to read graph number %d: %s \n",i,name[i]);

// Reading graph
//    if (!read_graph(&(n[i]),&(m[i]),name[i],&ver,&edges,&e,&w,&weight)) {
    if (!read_graph(&(n[i]),&(m[i]),name[i],&ver,&edges,&e,&w,&weight)) {
      printf("Problem reading graph %s \n",name[i]);
    }

  dist_e = (int *) malloc(2*sizeof(int)*(m[i]));
  for(jj=0;jj<m[i];jj++) { 
    dist_e[jj] = edges[jj];
  }

  left = (int *) malloc(sizeof(int)*n[i]);
  next = (int *) malloc(sizeof(int)*n[i]);
  order = (int *) malloc(sizeof(int)*n[i]);
  used = (int *) malloc(sizeof(int)*n[i]);

  if (left == NULL) {
    printf("Unable to allocate space for left-array \n");
    return(0);
  }

// List of best local suitor for each vertex in the weighted algorithm

  ws  = (float *) malloc(sizeof(float)*(*n));
  ws2 = (float *) malloc(sizeof(float)*(*n));
  if ((ws == NULL) || (ws2 == NULL)) {
    printf("Unable to allocate space for ws-arrays in read_graph \n");
    return(0);
  }

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
    printf("Sequential unweighted vertex based greedy matching took %f seconds\n",t_bs[i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      eseq(n[i],m[i],e,p);
      mt2 = omp_get_wtime(); 
      if ((t_bse[i] < 0) || (mt2-mt1 < t_bse[i]))
        t_bse[i] = mt2-mt1;
    }

    printf("Sequential unweighted edge based greedy matching took %f seconds\n",t_bse[i]);


    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      wseq(n[i],ver,edges,p,weight);
      mt2 = omp_get_wtime(); 
      if ((t_bgw[i] < 0) || (mt2-mt1 < t_bgw[i]))
        t_bgw[i] = mt2-mt1;
    }

    printf("Sequential local maximal weight matching took %f seconds\n",t_bgw[i]);

// Compute a random ordering of the vertices
    for(l=1;l<=n[i];l++) {
      order[l] = l;
    }
    srandom(1);
    for(l=1;l<n[i];l++) {
      long int x = random() % (long int) (n[i]-l+1);
      x++;
      int tmp = order[n[i]-l+1];
      order[n[i]-l+1] = order[x];
      order[x] = tmp;
    }

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      wseq2(n[i],ver,edges,p,weight,order);
      mt2 = omp_get_wtime(); 
      if ((t_rgw[i] < 0) || (mt2-mt1 < t_rgw[i]))
        t_rgw[i] = mt2-mt1;
    }

    printf("Random order local maximal weight matching took %f seconds\n",t_rgw[i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      wseq3(n[i],ver,edges,p,weight,order);
      mt2 = omp_get_wtime(); 
      if ((t_rgw[i] < 0) || (mt2-mt1 < t_rgw[i]))
        t_rgw[i] = mt2-mt1;
    }

    printf("Random order local maximal weight matching ** With Prefetching** took %f seconds\n",t_rgw[i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight(n[i],ver,edges,p,ws,weight);
      mt2 = omp_get_wtime(); 
      if ((t_bsw[i] < 0) || (mt2-mt1 < t_bsw[i]))
        t_bsw[i] = mt2-mt1;
    }

    printf("Sequential weighted matching took %f seconds\n",t_bsw[i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight3(n[i],ver,edges,p,ws,weight);
      mt2 = omp_get_wtime(); 
      if ((t_bsw[i] < 0) || (mt2-mt1 < t_bsw[i]))
        t_bsw[i] = mt2-mt1;
    }

    printf("New Sequential weighted matching took %f seconds\n",t_bsw[i]);


    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight4(n[i],ver,edges,p,ws,p2,ws2,weight,used);
      mt2 = omp_get_wtime(); 
      if ((t_dgw[i] < 0) || (mt2-mt1 < t_dgw[i]))
        t_dgw[i] = mt2-mt1;
    }

    printf("Two round sequential weighted matching took %f seconds\n",t_dgw[i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight5(n[i],ver,edges,p,ws,p2,ws2,weight,used);
      mt2 = omp_get_wtime(); 
      if ((t_5gw[i] < 0) || (mt2-mt1 < t_5gw[i]))
        t_5gw[i] = mt2-mt1;
    }

    printf("Two round sequential weighted matching with separate dynamic programming took %f seconds\n",t_5gw[i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight6(n[i],ver,edges,p,ws,p2,ws2,weight,used);
      mt2 = omp_get_wtime(); 
      if ((t_6gw[i] < 0) || (mt2-mt1 < t_6gw[i]))
        t_6gw[i] = mt2-mt1;
    }

    printf("Two round greedy weighted matching with separate dynamic programming took %f seconds\n",t_6gw[i]);

    wed *we;

    we = (wed *) malloc(m[i] * sizeof(wed));
    if (we == NULL)
      printf("Allocation for wed did not succeed!!! \n");

    for(l=0;l<m[i];l++) {
      we[l].x = e[l].x;
      we[l].y = e[l].y;
      we[l].w = w[l];
      if (we[l].w != w[l])
        printf("Data did not copy %d \n",l);
    }


    mt1 = omp_get_wtime(); 
    qsort(we,m[i],sizeof(wed),compedge);
    mt2 = omp_get_wtime(); 
    printf("Time to sort the edges is %f seconds \n",mt2-mt1);

    for(l=0;l<m[i]-1;l++) {
      if (we[l].w < we[l+1].w) {
        printf("Data is not sorted in decending order \n");
        printf("%d =  %f   < %d =  %f \n",l,we[l].w,l+1,we[l+1].w);
        exit(0);
      }
    }


    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      gpa(n[i],m[i],we,ver,edges,weight);
      mt2 = omp_get_wtime(); 
      if ((t_gpa[i] < 0) || (mt2-mt1 < t_gpa[i]))
        t_gpa[i] = mt2-mt1;
    }

    printf("GPA algorithm with dynamic programming took %f seconds\n",t_gpa[i]);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      path(n[i],ver,edges,weight);
      mt2 = omp_get_wtime(); 
      if ((t_pat[i] < 0) || (mt2-mt1 < t_pat[i]))
        t_pat[i] = mt2-mt1;
    }

    printf("Path growing algorithm with dynamic programming took %f seconds\n",t_pat[i]);


/*
    mt1 = omp_get_wtime(); 
    wsort(n[i],ver,edges,weight);
    mt2 = omp_get_wtime(); 
    printf("Sorting took %f seconds\n",mt2-mt1);

    for(l=0;l<n_runs;l++) {
      mt1 = omp_get_wtime(); 
      sweight1(n[i],ver,edges,p,ws,weight,next);
      mt2 = omp_get_wtime(); 
      if ((t_bsw[i] < 0) || (mt2-mt1 < t_bsw[i]))
        t_bsw[i] = mt2-mt1;
    }

    printf("Sequential sorted weighted matching took %f seconds\n",t_bsw[i]);
*/

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
      t1_edge[i][j] = -1;
      t_ver[i][j] = -1;
      t_lock[i][j] = -1;
      t_weight[i][j] = -1;
      t1_weight[i][j] = -1;
      t2_weight[i][j] = -1;

    

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

#pragma omp barrier

// Running optimized edge based verification code
/*
#pragma omp master
          {  mt1 = omp_get_wtime(); }
             edgec1(n[i],m[i],e,p,left);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t1_edge[i][j] < 0) || (mt2-mt1 < t1_edge[i][j]))
               t1_edge[i][j] = mt2-mt1;
          }
*/

#pragma omp barrier

// Running edge based locking code

#pragma omp master
          {  mt1 = omp_get_wtime(); }
             locking1(n[i],m[i],e,p,nlocks);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t_lock[i][j] < 0) || (mt2-mt1 < t_lock[i][j]))
               t_lock[i][j] = mt2-mt1;
          }

#pragma omp barrier

// Running weighted algorithm


#pragma omp master
          {  mt1 = omp_get_wtime(); }
             weighted(n[i],ver,edges,p,ws,left,nlocks,weight);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t_weight[i][j] < 0) || (mt2-mt1 < t_weight[i][j]))
               t_weight[i][j] = mt2-mt1;
          }
// Running second weighted algorithm
#pragma omp barrier

#pragma omp master
          {  mt1 = omp_get_wtime(); }
             weighted2(n[i],ver,edges,p,ws,nlocks,weight);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t1_weight[i][j] < 0) || (mt2-mt1 < t1_weight[i][j]))
               t1_weight[i][j] = mt2-mt1;
          }
// Running third weighted algorithm

#pragma omp barrier

         
// Move the data to a new data structure

  int start = my_id*n[i]/threads;
  int final = start + n[i]/threads;
  int ii,jj;

  for(ii=start;ii<=final;ii++) {                // Loop over vertices
    for(jj=ver[ii];jj<ver[ii+1];jj++) { // Loop over neighbors of vertex i
      dist_e[jj] = edges[jj];
    }
  } // loop over neighbors



#pragma omp barrier


#pragma omp master
          {  mt1 = omp_get_wtime(); }
             weighted3(n[i],m[i],ver,dist_e,p,ws,nlocks,weight,left);
#pragma omp master
          {  mt2 = omp_get_wtime(); 
             if ((t2_weight[i][j] < 0) || (mt2-mt1 < t2_weight[i][j]))
               t2_weight[i][j] = mt2-mt1;
          }


        } // End of k iterations over same configuration 
       
#pragma omp master
        { printf("Vertex code using %d threads took %f seconds, absolute speedup = %4.2f \n",conf[j],t_ver[i][j],t_bs[i]/t_ver[i][j]);
          printf("Vertex code using %d threads took %f seconds, relative speedup = %4.2f \n",conf[j],t_ver[i][j],t_ver[i][0]/t_ver[i][j]);
          printf("Edge code using %d threads took %f seconds, absolute speedup = %4.2f \n",conf[j],t_edge[i][j],t_bse[i]/t_edge[i][j]);
          // printf("** Optimized edge code using %d threads took %f seconds, absolute speedup = %4.2f \n",conf[j],t1_edge[i][j],t_bse[i]/t1_edge[i][j]);
          printf("Edge code using %d threads took %f seconds, relative speedup = %4.2f \n",conf[j],t_edge[i][j],t_edge[i][0]/t_edge[i][j]);
          printf("Locking code using %d threads took %f seconds, absolute speedup = %4.2f \n",conf[j],t_lock[i][j],t_bse[i]/t_lock[i][j]);
          printf("Locking code using %d threads took %f seconds, relative speedup = %4.2f \n",conf[j],t_lock[i][j],t_lock[i][0]/t_lock[i][j]);

          printf("Weighted code using %d threads took %f seconds, absolute speedup = %4.2f \n",conf[j],t_weight[i][j],t_bsw[i]/t_weight[i][j]);
          printf("Weighted code using %d threads took %f seconds, relative speedup = %4.2f \n",conf[j],t_weight[i][j],t_weight[i][0]/t_weight[i][j]);
          printf("New weighted code using %d threads took %f seconds, absolute speedup = %4.2f \n",conf[j],t1_weight[i][j],t_bsw[i]/t1_weight[i][j]);
          printf("New weighted code using %d threads took %f seconds, relative speedup = %4.2f \n",conf[j],t1_weight[i][j],t1_weight[i][0]/t1_weight[i][j]);
          printf("Dummy weighted code using %d threads took %f seconds, absolute speedup = %4.2f \n",conf[j],t2_weight[i][j],t_bs[i]/t2_weight[i][j]);
          printf("Dummy weighted code using %d threads took %f seconds, relative speedup = %4.2f \n",conf[j],t2_weight[i][j],t2_weight[i][0]/t2_weight[i][j]);

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
        (*w)[num_edges] = fabsf(v);

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

