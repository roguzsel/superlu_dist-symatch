// weighted algorithm
//
// Optimized main loop
// Two round sequential suitor matching followed by dynamic programming
//
// Note that this code does not return the actual matching, only its weight
// For this reason it is slightly faster than the other codes

void sweight4(int n,int *ver,int *edges,
              int *s,double *ws,		// Storage for first matching
              int *s2, double *ws2,		// Storage for second matching
              double *weight, int *used) {


  int i,j;
  int partner;     // Prospective candidate to match with   
  int x, done, next_vertex;
  int count,current;

// Start of matching algorithm

  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0

    s2[i] = 0;           // Set that no node is trying to match with i 
    ws2[i] = 0.0;        // Set the current weight of best suitor edge to 0

    used[i] = false;     // Set that node i has not been used in the dynamic programming part
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
  int nr_paths = 0;
  int nr_single = 0;

  int cum_length = 0;
  double cum_weight= 0.0;
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    if ((s2[i] == 0) && (s[i] != 0)) {  // Starting with a level 1 edge
      nr_paths++;
      int length = 1;
      int w_old = 0.0;
      int w_new = ws[i];
      double tmp;
      used[i] = true;
      current = s[i];
      
      while (true) {   // Standing at the endpoint of an edge from the level 1 matching
        used[current] = true;   // Set the current edge as used
        if (s2[current] == 0)   // If the path does not continue then move to next path
          break;

        if (ws2[current] + w_old > w_new) {  // Pick the best matching
          tmp = ws2[current] + w_old;
          w_old = w_new;
          w_new = tmp;
        }
        else {
          w_old = w_new;
        }

        current = s2[current];    // Move along the path using an edge from the level 2 matching
        used[current] = true;     // Mark new vertex as used
        length++;
        if (s[current] == 0)      // If the path ends, move to next path
          break;

        if (ws[current] + w_old > w_new) { // Pick best matching
          tmp = ws[current] + w_old;
          w_old = w_new;
          w_new = tmp;
        }
        else {
          w_old = w_new;
        }
        current = s[current]; // Move to the next vertex on the path
        length++;
      }
      cum_length += length;
      cum_weight += w_new;
    }
    else if ((s2[i] != 0) && (s[i] == 0)) {  // Starting with a level 2 edge
      nr_paths++;
      int length = 1;
      int w_old = 0.0;
      int w_new = ws2[i];
      double tmp;
      used[i] = true;
      current = s2[i];
      
      while (true) {   // Standing at the endpoint of an edge from the level 2 matching
        used[current] = true;   // Set the current edge as used
        if (s[current] == 0)   // If the path does not continue then move to next path
          break;

        if (ws[current] + w_old > w_new) {  // Pick the best matching
          tmp = ws[current] + w_old;
          w_old = w_new;
          w_new = tmp;
        }
        else {
          w_old = w_new;
        }

        current = s[current];    // Move along the path using an edge from the level 1 matching
        used[current] = true;     // Mark new vertex as used
        length++;
        if (s2[current] == 0)      // If the path ends, move to next path
          break;

        if (ws2[current] + w_old > w_new) { // Pick best matching
          tmp = ws2[current] + w_old;
          w_old = w_new;
          w_new = tmp;
        }
        else {
          w_old = w_new;
        }
        current = s2[current]; // Move to the next vertex on the path
        length++;
      }
      cum_length += length;
      cum_weight += w_new;
    }
    else if ((s2[i] == 0) && (s[i] == 0)) {  // Found single vertex
      used[i] = true;
      nr_single++;
    }
  }

// Now look for the cycles
  int nr_cycle = 0;
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    nr_cycle++;
    int length = 1;
    double tmp;

// Manually processing the first 5 vertices of the cycle, v1, v2, v3, v4, v5 (note v1 = v5 is a possibility)
// This is done since we must consider ignoring either (v1,v2) (using w1) or (v2,v3) (using w2)

    used[i] = true;		// Set v1 as used and move on, ignoring edge v1,v2 for now
    current = s[i];		// Move to v2
    used[current] = true;	// Set v2 as used
    current = s2[current];	// Move to v3
    used[current] = true;	// Set v3 as used

    double w1_old ,w1_new;  // Used when ignoring (v1,v2)
    double w2_old ,w2_new;  // Used when ignoring (v2,v3)

    if (ws[current] > ws2[current]) {  // Either pick (v3,v4) or (v2,v3)
      w1_new = ws[current];	// Pick (v3,v4);
    }
    else {			// Pick (v2,v3)
      w1_new = ws2[current];
    }
    w1_old = ws2[current];  

    current = s[current];	// Move to v4
    used[current] = true;	// Set v4 as used

// Ignoring (v2,v3). Pick best of (v3,v4) and (v4,v5)
    if (ws[current] > ws2[current]) {  // Either pick (v3,v4) or (v4,v5)
      w2_new = ws[current];	// Pick (v3,v4)
    }
    else {			// Pick (v4,v5)
      w2_new = ws2[current];
    }
    w2_old = ws[current];

    current = s2[current];		// Move to v5

    if (ws2[current] + w1_old > w1_new) {  // Pick the best matching
      tmp = ws2[current] + w1_old;      // Include (v4,v5)
      w1_old = w1_new;
      w1_new = tmp;
    }
    else {
      w1_old = w1_new;			// Ignore (v4,v5)
    }

// Ready to start loop
      
    while (!used[current]) {   // Process two edges (vi,v(i+1)) and (v(i+1),v(i+2)), at a time until we reach the start of the cycle
      used[current] = true;    // Set vi as used
      current = s[current];    // Move to v(i+1)
      used[current] = true;    // Set v(i+1) as used

// First do w1 values
      if (ws[current] + w1_old > w1_new) {  // Pick the best matching
        tmp = ws[current] + w1_old;
        w1_old = w1_new;
        w1_new = tmp;
      }
      else {
        w1_old = w1_new;
      }
// Then do w2 values
      if (ws[current] + w2_old > w2_new) {  // Pick the best matching
        tmp = ws[current] + w2_old;
        w2_old = w2_new;
        w2_new = tmp;
      }
      else {
        w2_old = w2_new;
      }

      current = s2[current];    // Move to v(i+2)

// First do w1 values
      if (ws2[current] + w1_old > w1_new) {  // Pick the best matching
        tmp = ws2[current] + w1_old;
        w1_old = w1_new;
        w1_new = tmp;
      }
      else {
        w1_old = w1_new;
      }
// Then do w2 values
      if (ws2[current] + w2_old > w2_new) {  // Pick the best matching
        tmp = ws2[current] + w2_old;
        w2_old = w2_new;
        w2_new = tmp;
      }
      else {
        w2_old = w2_new;
      }

    } // End while loop

// Now consider use of (v1,v2) for path that ignores (v2,v3), that is, the w2 values.
    if (ws[current] + w2_old > w2_new) {  // Pick the best matching
      w2_new = ws[current] + w2_old;
    }
    
// Finaly, add in the weight of the best solution
 
    if (w1_new > w2_new)
      cum_weight += w1_new;
    else
      cum_weight += w2_new;

  } // End of loop over vertices that looks for unprocessed cycles

  // printf("Found %d cycle nodes \n",nr_cycle);

  // printf("Found %d paths of total length %d, weight %f and %d single in two level algorithm \n",nr_paths/2,cum_length,cum_weight,nr_single);

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
