// weighted algorithm
//
// Optimized main loop
// Two round sequential matching followed by dynamic programming
// This time the dynamic programming is done in a separate routine 
// Note that the cycle DP is done in a slightly inefficient way. First with two separate calls to the DP routine
// to find which solution is the best and then with a separate call to set the actual values. This is fixed in
// later versions.

void sweight5(int n,int *ver,int *edges,int *s,double *ws,int *s2, double *ws2,double *weight, int *used,int *match) {


  int i,j;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;
  int count,current;
  int *path1;
  int *path2;
  int l1,l2;
  double *weight1;
  double *weight2;
  double cum_weight = 0.0;
  double dyn_prog();

// This could probably be moved to the main algorithm


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

// First round of matching

  i = 1;
  current = i;

  while (current != n+1) {
    double heaviest = ws[current];          // Find the heaviest candidate that does not have a better offer
    partner = s[current];                  // No point in trying a partner worse than the best suitor

    for(j=ver[current];j<ver[current+1];j++) { // Loop over neighbors of the current vertex 
      int y = edges[j];            // y is the neighbor of the current vertex
      if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(current,y) is the best so far, and if it is a better option for y
        heaviest = weight[j];      // Store the weight of the heaviest edge found so far
        partner = y;               // Store the name of the associated neighbor
      } // loop over neighbors
    }

    
    if (heaviest > 0) {            // True if there is a new partner
      if (s[partner] != 0) {       // True if partner already had a suitor
        next_vertex = s[partner];  // Pick up the old suitor and continue
      }
      else {
        i = i + 1;                 // Move to the next vertex
        next_vertex = i;
      }
      s[partner] = current;        // current is now the current suitor of s
      ws[partner] = heaviest;      // The weight of the edge (current,partner) 
    }
    else {
      i = i + 1;                   // Move to the next vertex
      next_vertex = i;
    }
    current = next_vertex;         // Update current vertex with the next one to work on
  } // loop over vertices


// Second round of matching

  i = 1;
  current = i;

  while (current != n+1) {
    double heaviest = ws2[current];          // Find the heaviest candidate that does not have a better offer
    partner = s2[current];                  // No point in trying a partner worse than the best suitor

    for(j=ver[current];j<ver[current+1];j++) { // Loop over neighbors of the current vertex 
      int y = edges[j];            // y is the neighbor of the current vertex
      // if (y == s[current]) continue;
      if ((weight[j]>heaviest) && (ws2[y]<weight[j]) && (y != s[current])) {// Check if w(current,y) is the best so far, and if it is a better option for y. It must also be different from the vertex selected in stage 1
        heaviest = weight[j];      // Store the weight of the heaviest edge found so far
        partner = y;               // Store the name of the associated neighbor
      } // loop over neighbors
    }

    
    if (heaviest > 0) {            // True if there is a new partner
      if (s2[partner] != 0) {       // True if partner already had a suitor
        next_vertex = s2[partner];  // Pick up the old suitor and continue
      }
      else {
        i = i + 1;                 // Move to the next vertex
        next_vertex = i;
      }
      s2[partner] = current;        // current is now the current suitor of s
      ws2[partner] = heaviest;      // The weight of the edge (current,partner) 
    }
    else {
      i = i + 1;                   // Move to the next vertex
      next_vertex = i;
    }
    current = next_vertex;         // Update current vertex with the next one to work on
  } // loop over vertices

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


// Finally look for any remaining free vertices and match these with their heaviest edge

  int k;
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


  free(weight1);
  free(weight2);
  free(path1);
  free(path2);


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

// Dynamic programming on a weighted path to find the heaviest matching
// The path contains l1 vertices, stored in path1[]
// and l1-1 edges, the weights of these are stored in match[].

double dyn_prog(int l1,int *path1,double *weight1,int *match) {

  double temp;
  double w_old = 0.0;
  double w_new = weight1[0];
  int i,k;
  int in[l1];
  int partner;

  in[0] = true;   // The 0'th edge is in the solution

  match[path1[0]] = 0;
  for(k=1;k<l1;k++) {
    match[path1[k]] = 0;
    if (weight1[k] + w_old > w_new) {
      temp = weight1[k] + w_old;
      w_old = w_new;
      w_new = temp;
      in[k] = true;   // The k'th edge is in the solution
    }
    else {
      w_old = w_new;
      in[k] = false;  // The k'th edge is out of the solution
    }
  }
  match[path1[l1]] = 0;

// Now match up the vertices that belong to the final solution
  k = l1-1;
  while (k >= 0) {
    if (in[k]) {
      match[path1[k]] = path1[k+1];
      match[path1[k+1]] = path1[k];
      k = k - 2;
    }
    else {
      k = k - 1;
    }
  }

  return w_new;  // Return cost of best matching

}
