// Two level path growing, each with DP. Then followed by one more round of DP
//

void pth1(int n,int *ver,int *edges,double *weight,int *match ) {

  int i,j;
  double heaviest;
  int best;
  int y;
  int length;
  int *path;
  double *pw;
  int *visited;
  double sdyn_prog();
  int partner;
  int *p1,*p2;
  double *wp1,*wp2;

// Start of matching algorithm

//  printf("Starting sequential matching algorithm \n");

  pw = (double *) malloc((n+1)* sizeof(double));
  path = (int *) malloc((n+1)* sizeof(int));
  visited = (int *) malloc((n+1)* sizeof(int));
  p1 = (int *) malloc((n+1)* sizeof(int));
  p2 = (int *) malloc((n+1)* sizeof(int));
  wp1 = (double *) malloc((n+1)* sizeof(double));
  wp2 = (double *) malloc((n+1)* sizeof(double));

  for(i=1;i<=n;i++) {   // Set that every node is unmatched
    visited[i] = false;
    match[i] = 0;
    p1[i] = 0;
    p2[i] = 0;
  }

  for(i=1;i<=n;i++) {                 // Loop over vertices
    if (visited[i])
      continue;
    int done = false;
    int next = i;
    length = -1;
    while (!done) {
      length++;
      path[length] = next;   
      visited[next] = true; 
      heaviest = -1.0;
      for(j=ver[next];j<ver[next+1];j++) { // Loop over neighbors of vertex next
        y = edges[j];
        if ((!visited[y]) && (weight[j]>heaviest)) {     // Check if vertex is unmatched and if it gives a new heaviest edge
          heaviest = weight[j];
          best = y;
        }
      }

      if (heaviest < 0.0) 
        done = true;
      else {
        next = best;
        pw[length] = heaviest;
      }
    } // End while

    if (length == 0)
      continue; // There was only one vertex and no edge

// Now do dynamic programming
    sdyn_prog(length,path,pw,p1);

  } // loop over vertices


  int k;

  for(i=1;i<=n;i++) {
    visited[i] = false;    // Reuse this in next round of the algorithm
    if (p1[i] == 0) {
      double best = -1.0;
      for(k=ver[i];k<ver[i+1];k++) {
        int y = edges[k];
        if ((p1[y] == 0) && (weight[k] > best)) {
          best = weight[k];
          partner = y;
        }
      }
      if (best > 0.0) {
        p1[i] = partner;
        p1[partner] = i;
        wp1[i] = best;
        wp1[partner] = best;
      }
    }
    else {	// Store the weight of the edge that vertex i is matched with, this is not efficient..
      for(k=ver[i];k<ver[i+1];k++) {
        int y = edges[k];
        if (p1[i] == y) {
          wp1[i] = weight[k];
          break;
        }
      }
      // printf("p[%d] = %d, but p[%d] = %d \n",i,p1[i],p1[i],p1[p1[i]]);
      // return;
    } // else
  } // for i

//  Now do second round of path growing, but ignoring the edges from the first matching

  for(i=1;i<=n;i++) {                 // Loop over vertices
    if (visited[i])
      continue;
    int done = false;
    int next = i;
    length = -1;
    while (!done) {
      length++;
      path[length] = next;   
      visited[next] = true; 
      heaviest = -1.0;
      for(j=ver[next];j<ver[next+1];j++) { // Loop over neighbors of vertex next
        y = edges[j];
// Check if y is not used in the previous matching, y is unmatched and if it gives a new heaviest edge
        if ((p1[next] != y) && (!visited[y]) && (weight[j]>heaviest)) {     
          heaviest = weight[j];
          best = y;
        }
      }

      if (heaviest < 0.0) 
        done = true;
      else {
        next = best;
        pw[length] = heaviest;
      }
    } // End while

    if (length == 0)
      continue; // There was only one vertex and no edge

// Now do dynamic programming
      sdyn_prog(length,path,pw,p2);

  } // loop over vertices

// Fill up empty slots, but ignore edges used in p1

  for(i=1;i<=n;i++) {
    visited[i] = false;    // Reuse this in next round of the algorithm
    if (p2[i] == 0) {
      double best = -1.0;
      for(k=ver[i];k<ver[i+1];k++) {
        int y = edges[k];
        if ((p1[i] != y) && (p2[y] == 0) && (weight[k] > best)) {
          best = weight[k];
          partner = y;
        }
      }
      if (best > 0.0) {
        p2[i] = partner;
        p2[partner] = i;
        wp2[i] = best;
        wp2[partner] = best;
      }
    }
    else {	// Store the weight of the edge that vertex i is matched with, this is not efficient..
      for(k=ver[i];k<ver[i+1];k++) {
        int y = edges[k];
        if (p2[i] == y) {
          wp2[i] = weight[k];
          break;
        }
      } // for k
        // printf("p[%d] = %d, but p[%d] = %d \n",i,p1[i],p1[i],p1[p1[i]]);
        // return;
    } // else
  } // for i


// ****
// Do final round of DP using the p1 and the p2 matchings, associated weights are stored in wp1 and wp2
// ****
  int l1;
  int current;

  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (visited[i]) continue;  // If this vertex has been used previously then skip it

    if ((p2[i] == 0) && (p1[i] != 0)) {  // Starting with a level 1 edge
      l1 = 0;                   // Keep track of the length of this path, starting from 0
      path[l1] = i;            // Store the first vertex
      pw[l1] = wp1[i];      // Store the weight of the first edge
      visited[i] = true;           // Mark the vertex as used
      current = p1[i];           // Move to the next vertex

      while (true) {   // Standing at the endpoint of an edge from the level 1 matching
        visited[current] = true;   // Set the current vertex as used
        l1++;                   // Increase the length of the path
        path[l1] = current;    // Store this vertex
        if (p2[current] == 0)   // If the path does not continue then move to next path
          break;
        pw[l1] = wp2[current];     // Store the weight of the next edge

        current = p2[current];    // Move along the path using an edge from the level 2 matching
        visited[current] = true;     // Mark new vertex as used
        l1++;
        path[l1] = current;      // Store the next vertex
        if (p1[current] == 0)      // If the path ends, move to next path
          break;
        pw[l1] = wp1[current];      // Store the weight of the next edge

        current = p1[current]; // Move to the next vertex on the path
      } // End while

      sdyn_prog(l1,path,pw,match);

    } // End if
    else if ((p2[i] != 0) && (p1[i] == 0)) {  // Starting with a level 2 edge
      l1 = 0;                   // Keep track of the length of this path, starting from 0
      path[l1] = i;            // Store the first vertex
      pw[l1] = wp2[i];     // Store the weight of the first edge
      visited[i] = true;           // Mark the vertex as used
      current = p2[i];          // Move to the next vertex

      while (true) {   // Standing at the endpoint of an edge from the level 2 matching
        visited[current] = true;   // Set the current vertex as used
        l1++;                   // Increase the length of the path
        path[l1] = current;    // Store this vertex
        if (p1[current] == 0)    // If the path does not continue then move to next path
          break;
        pw[l1] = wp1[current];      // Store the weight of the next edge

        current = p1[current];     // Move along the path using an edge from the level 2 matching
        visited[current] = true;     // Mark new vertex as used
        l1++;
        path[l1] = current;      // Store the next vertex
        if (p2[current] == 0)      // If the path ends, move to next path
          break;
        pw[l1] = wp2[current];     // Store the weight of the next edge

        current = p2[current]; // Move to the next vertex on the path
      } // End while

      sdyn_prog(l1,path,pw,match);

    }
    else if ((p2[i] == 0) && (p1[i] == 0)) {  // Found single vertex
      visited[i] = true;
    }
  }

// Now look for the cycles
  int nr_cycle = 0;
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (visited[i]) continue;  // If this vertex has been used previously then skip it

    nr_cycle++;
    current = i;
    l1 = 0;

// Ready to start loop

    while (!visited[current]) {   // Process two edges (vi,v(i+1)) and (v(i+1),v(i+2)), at a time until we reach the start of the cycle

// Put vertex v_i and edge (vi,v(i+1)) in the path

      visited[current] = true;             // Set vi as used
      path[l1] = current;              // Put vi in path
      pw[l1] = wp1[current];        // Put weight of (vi,v(i+1)) in the path

      current = p1[current];             // Move to v(i+1)

      l1++;                             // Increase number of vertices
      visited[current] = true;             // Set v(i+1) as used
      path[l1] = current;              // Put v(i+1) in path
      pw[l1] = wp2[current];       // Put weight of (v(i+1),v(i+2) in the path

      current = p2[current];            // Move to v(i+2)
      l1++;
    } // End while loop

    path[l1] = path[0];               // Put the first vertex last to complete the cycle


// Now do dynamic programming on the cycle 
    cycle_dyn_prog(l1,path,pw,match);

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

  free(pw);
  free(path);
  free(visited);
  free(p1);
  free(p2);
  free(wp1);
  free(wp2);
}
