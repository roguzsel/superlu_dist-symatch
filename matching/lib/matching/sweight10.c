// weighted algorithm
//
// Assuming sorted edge lists!

// Two round algorithm followed by DP.

// Optimized main loop
// Only uses one while-loop to control the entire algorithm.
// After processing one node, the algorithm either moves to the next (linearly) if there is no
// misplaced suitor-node or it picks up the suitor node for processing.

void sweight10(int n,int *ver,int *edges,int *s1,double *ws1,int *s2,double *ws2,double *weight,int *next,int *used,int *match) {

  int i,y,e;
  int next_vertex;
  int current;
  int *path1;
  double *weight1;
  int l1;

  path1 = (int *) malloc(n * sizeof(int));
  if (path1 == NULL) {
    printf("Unable to allocate memory for path1 in sweight10 \n");
    return;
  }
  weight1 = (double *) malloc(n * sizeof(double));
  if (weight1 == NULL) {
    printf("Unable to allocate memory for weight1 in sweight10 \n");
    return;
  }

// Start of matching algorithm


  for(i=0;i<=n;i++) {   
    s1[i] = 0;           // Set that no node is trying to match with i 
    ws1[i] = -1.0;       // Set the current weight of best suitor edge to 0

    s2[i] = 0;           // Set that no node is trying to match with i 
    ws2[i] = -1.0;       // Set the current weight of best suitor edge to 0

    next[i] = ver[i];   // Set where to start searching for the next candidate
    match[i] = 0;
    used[i] = false;
  }

// Round 1 of matching

  i = 1;
  current = i;


  while (current != n+1) {
    y = next[current];    // y points to the next possible neighbor (in the edge list of current) that current can match with
    e = edges[y];         // Get the id of the neighbor that we will try to match with

// Stop if there are no more edges or if we found a candidate
    while ((y < ver[current+1]) && ((weight[y] < ws1[e]) || ((weight[y] == ws1[e]) && (current < s1[e])))) {
      y++;              // Move to next position in edge list
      e = edges[y];     // Get id of neighbor
    }

    if (y < ver[current+1]) {      // Test if we found a possible partner

      next[current] = y+1;         // Set where to search from next time current needs a partner

      if (s1[e] != 0) {             // True if e already had a suitor
        next_vertex = s1[e];        // Pick up the old suitor
      }
      else {
        i = i + 1;                 // Move to the next vertex
        next_vertex = i;
      }
      s1[e] = current;              // current is now the suitor of e
      ws1[e] = weight[y];           // Store the weight of (i,e) as the weight given by the current suitor for e
    }
    else {
      i = i + 1;                   // Move to the next vertex
      next_vertex = i;
    }
    current = next_vertex;         // Update current vertex with the next one to work on
  } // loop over vertices


// Round 2 of matching, cannot use edges from round 1

  for(i=0;i<=n;i++) {   
    next[i] = ver[i];   // Set where to start searching for next candidate
    if ((s1[i] != 0) && (s1[s1[i]] != i)) {
      printf("%d is a suitor of %d, but %d is a suitor of %d \n",s1[i],i,s1[s1[i]],s1[i]);
      return;
    }
  }

  i = 1;
  current = i;

  while (current != n+1) {

    y = next[current];    // y points to the next possible neighbor (in the edge list of current) that current can match with
    e = edges[y];         // Get the id of the neighbor that we will try to match with

// Stop if there are no more edges or if we found a candidate
    while ((y < ver[current+1]) && ((weight[y] < ws2[e]) || (e == s1[current]) || ((weight[y] == ws2[e]) && (current < s2[e])))) {
      y++;              // Move to next position in edge list                       different from the one chosen in the first round
      e = edges[y];     // Get id of neighbor
    }

    if (y < ver[current+1]) {      // Test if we found a possible partner

      next[current] = y+1;         // Set where to search from next time current needs a partner

      if (s2[e] != 0) {             // True if e already had a suitor
        next_vertex = s2[e];        // Pick up the old suitor
      }
      else {
        i = i + 1;                 // Move to the next vertex
        next_vertex = i;
      }
      s2[e] = current;              // current is now the suitor of e
      ws2[e] = weight[y];           // Store the weight of (i,e) as the weight given by the current suitor for e
    }
    else {
      i = i + 1;                   // Move to the next vertex
      next_vertex = i;
    }
    current = next_vertex;         // Update current vertex with the next one to work on
  } // loop over vertices

// Starting dynamic programming

  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if ((s2[i] != 0) && (s2[s2[i]] != i)) {
      printf("%d is matched to %d, but %d is matched to %d \n",i,s2[i],s2[i],s2[s2[i]]);
      return;
    }
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    if ((s2[i] == 0) && (s1[i] != 0)) {  // Starting with a level 1 edge
      l1 = 0;                   // Keep track of the length of this path, starting from 0
      path1[l1] = i;            // Store the first vertex
      weight1[l1] = ws1[i];      // Store the weight of the first edge
      used[i] = true;           // Mark the vertex as used
      current = s1[i];           // Move to the next vertex

      while (true) {   // Standing at the endpoint of an edge from the level 1 matching
        used[current] = true;   // Set the current vertex as used
        l1++;                   // Increase the length of the path
        path1[l1] = current;    // Store this vertex
        if (s2[current] == 0)   // If the path does not continue then move to next path
          break;
        weight1[l1] = ws2[current];     // Store the weight of the next edge

        current = s2[current];    // Move along the path using an edge from the level 2 matching
        used[current] = true;     // Mark new vertex as used
        l1++;
        path1[l1] = current;      // Store the next vertex
        if (s1[current] == 0)      // If the path ends, move to next path
          break;
        weight1[l1] = ws1[current];      // Store the weight of the next edge

        current = s1[current]; // Move to the next vertex on the path
      } // End while

      dyn_prog(l1,path1,weight1,match);

    } // End if
    else if ((s2[i] != 0) && (s1[i] == 0)) {  // Starting with a level 2 edge
      l1 = 0;                   // Keep track of the length of this path, starting from 0
      path1[l1] = i;            // Store the first vertex
      weight1[l1] = ws2[i];     // Store the weight of the first edge
      used[i] = true;           // Mark the vertex as used
      current = s2[i];          // Move to the next vertex

      while (true) {   // Standing at the endpoint of an edge from the level 2 matching
        used[current] = true;   // Set the current vertex as used
        l1++;                   // Increase the length of the path
        path1[l1] = current;    // Store this vertex
        if (s1[current] == 0)    // If the path does not continue then move to next path
          break;
        weight1[l1] = ws1[current];      // Store the weight of the next edge

        current = s1[current];     // Move along the path using an edge from the level 2 matching
        used[current] = true;     // Mark new vertex as used
        l1++;
        path1[l1] = current;      // Store the next vertex
        if (s2[current] == 0)      // If the path ends, move to next path
          break;
        weight1[l1] = ws2[current];     // Store the weight of the next edge

        current = s2[current]; // Move to the next vertex on the path
      } // End while

      dyn_prog(l1,path1,weight1,match);

    }
    else if ((s2[i] == 0) && (s1[i] == 0)) {  // Found single vertex
      used[i] = true;
    }
  }


// Now look for the cycles
  int nr_cycle = 0;
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    nr_cycle++;
    current = i;
    l1 = 0;

// Ready to start loop

    while (!used[current]) {   // Process two edges (vi,v(i+1)) and (v(i+1),v(i+2)), at a time until we reach the start of the cycle

// Put vertex v_i and edge (vi,v(i+1)) in the path

      used[current] = true;             // Set vi as used
      path1[l1] = current;              // Put vi in path1
      weight1[l1] = ws1[current];        // Put weight of (vi,v(i+1)) in the path

      current = s1[current];             // Move to v(i+1)

      l1++;                             // Increase number of vertices
      used[current] = true;             // Set v(i+1) as used
      path1[l1] = current;              // Put v(i+1) in path1
      weight1[l1] = ws2[current];       // Put weight of (v(i+1),v(i+2) in the path

      current = s2[current];            // Move to v(i+2)
      l1++;
    } // End while loop

    path1[l1] = path1[0];               // Put the first vertex last to complete the cycle


// Now do dynamic programming on the cycle 
    cycle_dyn_prog(l1,path1,weight1,match);

  } // End of loop over vertices that looks for unprocessed cycles

  int k,partner;
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
  free(weight1);

}
