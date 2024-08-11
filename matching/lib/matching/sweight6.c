// weighted algorithm
//
// 
// Two round matching followed by dynamic programming
// Using simple locally greedy matching
// This version is using the inefficient cycle-DP. Don't know how much this hurts...

void sweight6(int n,int *ver,int *edges,int *s,double *ws,int *s2, double *ws2,double *weight, int *used,int *match) {


  int i,j,k;
  int partner;     // Prospective candidate to match with   
  int x, y, done, next_vertex;
  int count,current;
  int *path1;
  int *path2;
  int l1,l2;
  double *weight1;
  double *weight2;
  double cum_weight = 0.0;
  double dyn_prog();
  double heaviest;

  path1 = (int *) malloc(n * sizeof(int));
  path2 = (int *) malloc(n * sizeof(int));
  weight1 = (double *) malloc(n * sizeof(double));
  weight2 = (double *) malloc(n * sizeof(double));

// Start of matching algorithm

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0

    s2[i] = 0;           // Set that no node is trying to match with i 
    ws2[i] = 0.0;        // Set the current weight of best suitor edge to 0

    used[i] = false;     // Set that node i has not been used in the dynamic programming part
    match[i] = 0;
  }
//  count = 0;
  int num_match = 0;

// First round of matching

  for(i=1;i<=n;i++) {                 // Loop over vertices
    if (s[i] == 0) {                // Check if this vertex has been matched
      heaviest = 0.0;
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        y = edges[j];
        if ((s[y] == 0) && (weight[j]>heaviest)) {     // Check if vertex is unmatched and if it gives a new heaviest edge
          heaviest = weight[j];
          partner = y;
        }
      } // Loop over neighbors of vertex i

      if (heaviest > 0.0) {
        num_match++;
        s[i] = partner;           // Match i and partner
        s[partner] = i; 	     // Match both ways
        ws[i] = heaviest;      // Store the weight of the matched edge
        ws[partner] = heaviest;
      }

    }
  } // loop over vertices

  // printf("Number matched in first matching is %d \n",num_match);
  num_match = 0;

// Now do second matching, but ignoring edges used in the first one

  for(i=1;i<=n;i++) {                 // Loop over vertices
    if (s2[i] == 0) {                // Check if this vertex has been matched
      heaviest = 0.0;
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        y = edges[j];
        if ((s2[y] == 0) && (y != s[i]) && (weight[j]>heaviest)) {     // Check if vertex is unmatched, has not been used in the previous matching, and if it gives a new heaviest edge
          heaviest = weight[j];
          partner = y;
        }
      } // Loop over neighbors of vertex i

      if (heaviest > 0.0) {
        num_match++;
        s2[i] = partner;         // Match i and partner
        s2[partner] = i; 	     // Match both ways
        ws2[i] = heaviest;    // Store the weight of the matched edge
        ws2[partner] = heaviest;
      }

    }
  } // loop over vertices

  // printf("Number matched in second matching is %d \n",num_match);
// Ready to do find the paths and cycles and then do dynamic programming.

// Starting dynamic programming

  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    if ((s2[i] == 0) && (s[i] != 0)) {  // Starting with a level 1 edge
      l1 = 0;  			// Keep track of the length of this path, starting from 0
      path1[l1] = i;  		// Store the first vertex
      weight1[l1] = ws[i];	// Store the weight of the first edge
      used[i] = true;		// Mark the vertex as used
      current = s[i];   	// Move to the next vertex
      
      while (true) {   // Standing at the endpoint of an edge from the level 1 matching
        used[current] = true;   // Set the current vertex as used
        l1++;			// Increase the length of the path
        path1[l1] = current;  	// Store this vertex
        if (s2[current] == 0)   // If the path does not continue then move to next path
          break;
        weight1[l1] = ws2[current];	// Store the weight of the next edge

        current = s2[current];    // Move along the path using an edge from the level 2 matching
        used[current] = true;     // Mark new vertex as used
        l1++;
        path1[l1] = current;      // Store the next vertex
        if (s[current] == 0)      // If the path ends, move to next path
          break;
        weight1[l1] = ws[current];	// Store the weight of the next edge

        current = s[current]; // Move to the next vertex on the path
      } // End while

      dyn_prog(l1,path1,weight1,match);

    } // End if
    else if ((s2[i] != 0) && (s[i] == 0)) {  // Starting with a level 2 edge
      l1 = 0;  			// Keep track of the length of this path, starting from 0
      path1[l1] = i;  		// Store the first vertex
      weight1[l1] = ws2[i];	// Store the weight of the first edge
      used[i] = true;		// Mark the vertex as used
      current = s2[i];   	// Move to the next vertex
      
      while (true) {   // Standing at the endpoint of an edge from the level 2 matching
        used[current] = true;   // Set the current vertex as used
        l1++;			// Increase the length of the path
        path1[l1] = current;  	// Store this vertex
        if (s[current] == 0)    // If the path does not continue then move to next path
          break;
        weight1[l1] = ws[current];	// Store the weight of the next edge

        current = s[current];     // Move along the path using an edge from the level 2 matching
        used[current] = true;     // Mark new vertex as used
        l1++;
        path1[l1] = current;      // Store the next vertex
        if (s2[current] == 0)      // If the path ends, move to next path
          break;
        weight1[l1] = ws2[current];	// Store the weight of the next edge

        current = s2[current]; // Move to the next vertex on the path
      } // End while

      dyn_prog(l1,path1,weight1,match);

    }
    else if ((s2[i] == 0) && (s[i] == 0)) {  // Found single vertex
      used[i] = true;
    }
  }

// Now look for the cycles
  int nr_cycle = 0;
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    nr_cycle++;
    double tmp;

// Manually processing the first 5 vertices of the cycle, v1, v2, v3, v4, v5 (note v1 = v5 is a possibility)
// This is done since we must consider ignoring either (v1,v2) (using w1) or (v2,v3) (using w2)

    used[i] = true;		// Set v1 as used and move on

    current = s[i];		// Move to v2
    used[current] = true;	// Set v2 as used
    path1[0] = current;		// Put v2 in path1
    weight1[0] = ws2[current];

    current = s2[current];	// Move to v3
    used[current] = true;	// Set v3 as used
    path1[1] = current;		// Put v3 in path1
    weight1[1] = ws[current];
    path2[0] = current;		// Put v3 in path2
    weight2[0] = ws[current];

    current = s[current];	// Move to v4
    used[current] = true;	// Set v4 as used
    path1[2] = current;		// Put v4 in path1
    weight1[2] = ws2[current];
    path2[1] = current;		// Put v4 in path2
    weight2[1] = ws2[current];

    current = s2[current];	// Move to v5

    path1[3] = current;		// Put v5 in path1
    weight1[3] = ws[current];
    path2[2] = current;		// Put v5 in path2
    weight2[2] = ws[current];

    l1 = 3;
    l2 = 2;

// Ready to start loop
      
    while (!used[current]) {   // Process two edges (vi,v(i+1)) and (v(i+1),v(i+2)), at a time until we reach the start of the cycle
                               // Note that vi is already in both path1 and path2 (but not the weight of (vi,v(i+1))
      used[current] = true;    // Set vi as used
      weight1[l1] = ws[current];
      weight2[l2] = ws[current];

      current = s[current];    		// Move to v(i+1)
      used[current] = true;    		// Set v(i+1) as used
      l1++;
      path1[l1] = current;		// Put v(i+1) in path1
      weight1[l1] = ws2[current];
      l2++;
      path2[l2] = current;		// Put v(i+1) in path2
      weight2[l2] = ws2[current];

      current = s2[current]; 		// Move to v(i+2)
      l1++;
      path1[l1] = current;		// Put v(i+2) in path1
      l2++;
      path2[l2] = current;		// Put v(i+2) in path2
    } // End while loop

// Now add in (v1,v2) to path2
    weight2[l2] = ws[current];
    l2++;
    path2[l2] = s[current];		// Put v(i+1) in path2

// Now check which is better of path1 and path2. Redo the computation so that the best is left
    if (dyn_prog(l1,path1,weight1,match) > dyn_prog(l2,path2,weight2,match))
      dyn_prog(l1,path1,weight1,match);
    else dyn_prog(l2,path2,weight2,match);
 
  } // End of loop over vertices that looks for unprocessed cycles


// Check if there are any unmatched vertices as a result of the DP

  for(i=1;i<=n;i++) {
    if (match[i] == 0) {
      double best = -1.0;
      for(k=ver[i];k<ver[i+1];k++) {
        int y = edges[k];
        if ((match[y] == 0) && (weight[k] > best)) {
          best = weight[k];
          partner = y;
        }
      }
      if (best > 0.0) {
        match[i] = partner;
        match[partner] = i;
      }
    }
  }

  free(path1);
  free(path2);
  free(weight1);
  free(weight2);

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

