#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <time.h>


#define FALSE 0
#define TRUE  1

#define max_graph 28158+1

int main(int argc, char *argv[]) {

  int *weight[max_graph];
  int matched[max_graph];
  int i,j;
  int n;
  float mt1,mt2,tt;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;

  int s[max_graph];
  int ws[max_graph];

  n = max_graph-1;

  for(i=1;i<=n;i++) {
    matched[i] = FALSE;
    weight[i] = (int *) malloc(sizeof(int)*(n+1));
  }

  for(i=1;i<=n;i++) {
    weight[i][i] = -1;
    for(j=i+1;j<=n;j++) {
      weight[i][j] = rand();
      weight[j][i] = weight[i][j];
    }
  }


  mt1 = omp_get_wtime();

// Start of matching algorithm

  for(i=0;i<=n;i++) {
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0;          // Set the current weight of best suitor edge to 0
  }

  for(x=1;x<=n;x++) {                // Loop over vertices
    int done = FALSE;
    i = x;
    while (!done) {
      int heaviest = 0;
      for(j=1;j<=n;j++) { // Loop over neighbors of vertex i
        if ((weight[i][j]>heaviest) && (ws[j]<weight[i][j])) {// Check if w(i,y) is the best so far, and if it is a better option for j
          heaviest = weight[i][j];      // Store the best so far
          partner = j;
        }
      } // loop over neighbors

      done = TRUE;
      if (heaviest > 0) {
        if (s[partner] != 0) {  // True if partner already had a suitor
          next_vertex = s[partner];  // Pick up the old suitor and continue
          done = FALSE;
        }
        s[partner] = i;          // i is now the current suitor of s
        ws[partner] = heaviest;  // The weight of the edge (i,partner) 
      }
      if (!done) {
        i = next_vertex;
      }
    } // while not done
  } // loop over vertices

  mt2 = omp_get_wtime();
  tt = mt2-mt1;

  printf("Sequential greedy matching took %f seconds\n",tt);
  printf("Node 1 is matched with %d and is of weight %d \n",s[1],weight[1][s[1]]);

  for(x=1;x<=n;x++) {                // Loop over vertices
    if ((s[x] < 1) || (s[x] > n))
      printf("Error, s[%d] = %d  \n",x,s[x]);
    if (s[s[x]] != x) 
      printf("%d wants %d, but %d wants %d \n",x,s[x],s[x],s[s[x]]);
  } 

}
