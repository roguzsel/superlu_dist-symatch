// weighted algorithm
//
// This code uses a for-loop as the outer loop, with a while-loop to controll 
// the processing of suitor-nodes that need to be updated.
// This is the type of code that is used in the parallel version

void sweight(int n,int *ver,int *edges,int *s,double *ws,double *weight,int *init) {

	/* printf("n = %d\n", n); */
	/* int v; */
	/* printf("adj-list:\n"); */
	/* for (v = 1; v <= n; ++v) */
	/* { */
	/* 	printf("%3d : ", v); */
	/* 	for (int eptr = ver[v]; eptr < ver[v+1]; ++eptr) */
	/* 	{ */
	/* 		int e = edges[eptr]; */
	/* 		double w = weight[eptr]; */
	/* 		printf("(%5d, %.6lf) ", e, w); */
	/* 	} */
	/* 	printf("\n"); */
	/* } */
	

  int i,j;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;
  int count,ecount;

// Start of matching algorithm


/*
  int *init;
  free(init);
  init  = (int *) malloc((n+1)*sizeof(int));
  if (init == NULL) {
    printf("Unable to allocate memory for init in sweight() \n");
    return;
  }
*/

  // printf("AAA\n"); fflush(stdout);

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
    init[i] = 0;
  }

//   Initialize matching
  double mt1 = omp_get_wtime();


  for(x=1;x<=n;x++) {                // Loop over vertices
      double heaviest = ws[x];       // Find the heaviest candidate that does not have a better offer
      partner = s[x];                // No point in trying a partner worse than the best suitor
      for(j=ver[x];j<ver[x+1];j++) { // Loop over neighbors of vertex x
        int y = edges[j];            // y is the current neighbor of x

        if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y
          heaviest = weight[j];      // Store the weight of the heaviest edge found so far
          partner = y;               // Store the name of the associated neighbor
        }
      } // loop over neighbors

      if (heaviest > 0.0) {          // True if there is a partner
        if (s[partner] != 0) {
          init[s[partner]] = false;
        }
        s[partner] = x;              // x is now the current suitor of s
        ws[partner] = heaviest;      // The weight of the edge (x,partner) 
        init[x] = true;

        // printf("%d is matched to %d for a cost of %f \n",x,partner,heaviest);
      }


  }
 //  printf("Done initializing \n");
//

  double mt2 = omp_get_wtime();

//  printf("Time to initialize is %f \n",mt2-mt1);

  count = 0;
  ecount = 0;

  for(x=1;x<=n;x++) {                // Loop over vertices
   //  printf("Considering %d \n",x);
    if (init[x])
      continue;

   // if (x==1691)
   //   printf("%d initially wanted %d, now the suitor of %d is %d \n",x,init[x],init[x],s[init[x]]);

    int done = false;
    i = x;                           // Looking for a matching partner for node i
    while (!done) {

      count++;
      double heaviest = ws[i];          // Find the heaviest candidate that does not have a better offer
      partner = s[i];                  // No point in trying a partner worse than the best suitor
      
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        ecount++;
        int y = edges[j];            // y is the current neighbor of i

        if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y

/*
// This type of check is only needed in the parallel code

        if ((weight[j] < heaviest) || (weight[j] < ws[y]))
          continue;

        if ((weight[j] == heaviest) && (y < partner))
          continue;

        if ((weight[j] == ws[y]) && (i < s[y]))
          continue;
*/
          heaviest = weight[j];      // Store the weight of the heaviest edge found so far
          partner = y;               // Store the name of the associated neighbor
        }
      } // loop over neighbors

      done = true;                   // Check if we found someone to match with

      if (heaviest > 0) {            // True if there is a partner
        if (s[partner] != 0)  {       // True if partner already had a suitor
          next_vertex = s[partner];  // Pick up the old suitor and continue
      //     if (x==1691)
      //     printf("Picking up %d, increase in weight is %f \n",next_vertex,heaviest - ws[partner]);
          done = false;
        }
        s[partner] = i;              // i is now the current suitor of s
        ws[partner] = heaviest;      // The weight of the edge (i,partner) 
     //   if (x==1691)
     //    printf("%d is matched to %d for a cost of %f \n",i,partner,heaviest);
      }
      if (!done) {
        i = next_vertex;
      }
    } // while not done
  // printf("Vertex %d, Chain length = %d \n",x,count);
  } // loop over vertices



//  printf("Number of vertices considered %d, number of edges %d \n",count,ecount);
/*
   printf("Average number of iterations %f \n",(double)count/(double)n);

  int problems = 0;
  for(i=1;i<=n;i++) {
    if ((s[i] != 0) && (s[s[i]] != i))
      problems++;
  }
  if (problems > 0)
    printf("Have %d nodes with problems\n",problems);
*/
}
