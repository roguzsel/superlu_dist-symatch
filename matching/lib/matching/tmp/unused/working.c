#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "frame1.h"


int main(int argc, char *argv[]) {

  int i,j;

  FILE *rf;
  char *name;       // Name of current graph
  int n_graphs;     // Number of graphs to run
  int n_runs;       // Number of runs for each configuration
  int n_conf;       // Number of thread configurations
  int conf[100];    // Number of threads in each configuration

  double mt1,mt2;

  if (argc != 2) {
    printf("Give data file name as parameter\n");
    return;
  }
  printf("starting reading from %s \n",argv[1]);
// Opening file containing file names for graphs

  rf = fopen(argv[1],"r");
  if (rf == NULL) {
    printf("Cannot open data file \n");
    return;
  }

// Get number of graphs to read
  fscanf(rf,"%d",&n_graphs);
// Get number of runs per configuration
  fscanf(rf,"%d",&n_runs);
// Get number of thread configurations
  fscanf(rf,"%d",&n_conf);
// Get the different configurations

  for(i=0;i<n_conf;i++) {
    fscanf(rf,"%d",&(conf[i]));
  }
 
// Read name of graph
    fscanf(rf,"%s",name);

  fclose(rf);
  printf("Normal exit\n");
}

