#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int main(int argc, char *argv[]) {

//  int read_graph();
//  void discontinue();
  int i,j;

  FILE *rf;
  FILE *wf;
  char name[100];       // Name of current graph
  char *tc;         // Name of current graph without .mtx
  int n_graphs;     // Number of graphs to run
  int n_runs;       // Number of runs for each configuration
  int n_conf;       // Number of thread configurations
  int conf[100];    // Number of threads in each configuration



// Opening file containing file names for graphs
  rf = fopen("data_files","r");
/*
  n_graphs = 1;
  n_runs = 3;
  n_conf = 6;
*/

// Get the different configurations
  for(i=0;i<6;i++) {
    fscanf(rf,"%d",&(conf[i]));
  }

  for(i=0;i<6;i++) 
    printf("%d ",conf[i]);
  printf("\n");

  printf("starting reading of file name \n");
 
// Read name of graph
  fscanf(rf,"%s",name);

  printf("done reading of file name \n");

  fclose(rf);

  printf("Normal exit\n");
}

