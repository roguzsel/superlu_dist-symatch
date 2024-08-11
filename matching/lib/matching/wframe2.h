#ifndef _WFRAME2_H_
#define _WFRAME2_H_

#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* typedef struct { int x,y; } edge; */

/* typedef struct { int x,y; */
/*                  double w; } wed; */

/* typedef struct { int x; */
/*                  double w; } neig; */

/* #define true  1 */
/* #define false 0 */
/* #define TRUE  1 */
/* #define FALSE 0 */

/* /\* #define max_n  68000000 *\/ */
/* /\* #define max_m 761000000 *\/ */
/* #define max_n  68000 */
/* #define max_m 761000 */

/* #define max_experiment 200   // Maximum different algorithms that can be tested */
/* #define max_conf 20          // Maximum different thread configurations that can be tested */

/* #define max_graphs 20   // Maximum number of graphs in input file */

/* #define max_threads 500 */

int done = false;    // Global variable, used for signalling if some thread is not done

// Variables related to the setup specified in the input file

int n_graphs;        // Number of input graphs to run
int n_runs;          // Number of runs on each graph for each configuration
int n_conf;          // Number of thread configurations
int *conf;           // Number of threads in each configuration
char **name;         // Table with the names of the graphs 

FILE *wf;            // File pointer to output file

int *m; 	// Number of edges
int *n; 	// Number of vertices

double **timer;    // Timers for different runs
double **cost;    // Cost of solution for different runs
double ***p_timer; // Timers for parallel runs
double ***p_cost;  // Cost for parallel runs


#endif
