// weighted algorithm
//
// This code uses a for-loop as the outer loop, with a while-loop to controll 
// the processing of suitor-nodes that need to be updated.
// This is the type of code that is used in the parallel version.
//
// The algorithm continously removes any edge that cannot be used in the current matching. Thus
// it destroys the structure of  "edges" and "weight". It is assumed that nn is a copy of "ver".

void sweight11(int n,int *ver,int *edges,int *s,double *ws,double *weight,int *nn) {


  int i,j;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;
  int count,ecount;

// Start of matching algorithm
  int *init = (int *) malloc((n+1)*sizeof(int));

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0
    init[i] = 0;
  }


//   Initialize matching
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


  count = 0;
  ecount = 0;

  for(x=1;x<=n;x++) {                // Loop over vertices
    if (init[x])
      continue;

    int done = false;
    i = x;                           // Looking for a matching partner for node i
    while (!done) {

      count++;
      double heaviest = ws[i];         // Find the heaviest candidate that does not have a better offer
      partner = s[i];                  // No point in trying a partner worse than the best suitor
      

      for(j=ver[i];j<nn[i+1];j++)  { // Loop over neighbors of vertex i
        int y = edges[j];            // y is the current neighbor of i
        ecount++;
        if ((weight[j] < ws[i]) || (ws[y] >= weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y
/*
          int t1 = edges[j];
          double w1 = weight[j];
*/
          edges[j] = edges[nn[i+1]-1]; 
          weight[j] = weight[nn[i+1]-1]; 
/*
          edges[nn[i+1]-1] = t1; 
          weight[nn[i+1]-1] = w1; 
*/

          nn[i+1]--;
          j--;
        } else if (weight[j]>heaviest) {// Check if w(i,y) is the best so far

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


   // printf("Number of vertices considered %d, number of edges %d \n",count,ecount);
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
