// The GPA algorithm by Maue and Sanders

gpa(int n,int m,wed *we,int *ver,int *edges,float* weight,int *match) {

  int *parity;
  int *other;
  int *right;
  int *left;
  float *wr;
  float *wl;
  int v,w;
  int *path1;
  int *path2;
  float *weight1;
  float *weight2;
  float dyn_prog();
  int i;
  int current, next;
  int tmp;
  int length1,length2;
  int partner;

  parity = (int *) malloc(n * sizeof(int));
  other  = (int *) malloc(n * sizeof(int));
  right  = (int *) malloc(n * sizeof(int));
  left   = (int *) malloc(n * sizeof(int));
  wr     = (float *) malloc(n * sizeof(float));
  wl     = (float *) malloc(n * sizeof(float));

  path1   = (int *) malloc(n * sizeof(int));
  path2   = (int *) malloc(n * sizeof(int));
  weight1 = (float *) malloc(n * sizeof(float));
  weight2 = (float *) malloc(n * sizeof(float));

// Iterate through the edges (in order). For each edge check if it completes an even length cycle or
// if it connects two paths.

// First initialize variables for each vertex
  for(i=1;i<=n;i++) {
    parity[i] = true;
    other[i] = i;
    right[i] = 0;     
    left[i] = 0;     
    match[i] = 0;     
  }

// Then iterate through the edges

  for(i=1;i<=m;i++) {
    v = we[i].x;     
    w = we[i].y;     

// If either of v and w is not a path endpoint then skip this edge
    if (((right[v] != 0) && (left[v] != 0)) || ((right[w] != 0) && (left[w] != 0)))
      continue;

// Check if (v,w) completes a cycle
    if (other[v] == w) {  // we have a cycle

      if (parity[v])  // If we create an odd cycle then skip to next edge
        continue;

      parity[v] = true; // Set that this cycle is even
      parity[w] = true;

// Create the cycle
      if (right[v] == 0) {  // First link v to w, must find which pointer to use
        right[v] = w;
        wr[v] = we[i].w;  // Store the weight of the right edge of v
      }
      else {
        left[v] = w;
        wl[v] = we[i].w;  // Store the weight of the left edge of v
      }

      if (right[w] == 0) {  // Now connect w to v
        right[w] = v;
        wr[w] = we[i].w;
      }
      else {
        left[w] = v;
        wl[w] = we[i].w;
      }
      other[v] = -1;  // Set that this is a cycle
      continue;
    }  // End treatment of cycle

// (v,w) does not complete a cycle, v and w must be endpoints, ok to merge lists
        
    parity[other[v]] = parity[v]^parity[w];  // Performing logical XOR on the parity
    parity[other[w]] = parity[other[v]];     // Storing the new parity with the new endpoints 

    int tmp = other[v];
    other[tmp] = other[w];  // Update the remaining endpoints
    other[other[w]] = tmp;

// Link the lists, first v to w and then w to v
    if (right[v] == 0) {
      right[v] = w;
      wr[v] = we[i].w;
    }
    else {
      left[v] = w;
      wl[v] = we[i].w;
    }
    if (right[w] == 0) {
      right[w] = v;
      wr[w] = we[i].w;
    }
    else {
      left[w] = v;
      wl[w] = we[i].w;
    }

  } // End loop over m edges


// Now pick out the paths and cycles and do dynamic programming on them
  for(i=1;i<=n;i++) {
    if (other[i] == -1) {// This is a cycle
     
// Processing the first 5 vertices v1,v2,v3,v4,v5 manually, initially i = v1

      if (right[i] == 0)
        current = left[i];
      else
        current = right[i];

// Now current = v2

      path1[0] = current;

      if (right[current] == i) {
        next = left[current];
        weight1[0] = wl[current];  // Picking up weight of edge (v2,v3)
      }
      else {
        next = right[current];
        weight1[0] = wr[current];
      }

      tmp = current;
      current = next;

// Now current = v3 while tmp = v2

      path1[1] = current;
      path2[0] = current;

      if (right[current] == tmp) { // Moving next forward to v4
        next = left[current];
        weight1[1] = wl[current];  // Picking up weight of edge (v3,v4)
        weight2[0] = wl[current];  // Picking up weight of edge (v3,v4)
      }
      else {
        next = right[current];
        weight1[1] = wr[current];
        weight1[0] = wr[current];
      }

      tmp = current;
      current = next;

// Now current = v4 while tmp = v3

      path1[2] = current;
      path2[1] = current;

      if (right[current] == tmp) { // Moving next forward to v5
        next = left[current];
        weight1[2] = wl[current];  // Picking up weight of edge (v4,v5)
        weight2[1] = wl[current];  // Picking up weight of edge (v4,v5)
      }
      else {
        next = right[current];
        weight1[2] = wr[current];
        weight2[1] = wr[current];
      }

      tmp = current;
      current = next;

// Now current = v5 while tmp = v4

      length1 = 2;
      length2 = 1;

// Pick up rest of cycle
      while (current != i) {
        length1++;
        length2++;
        path1[length1] = current;
        path2[length2] = current;
        if (right[current] == tmp) { // Moving next forward to v5
          next = left[current];
          weight1[length1] = wl[current];  // Picking up weight of edge (v4,v5)
          weight2[length2] = wl[current];  // Picking up weight of edge (v4,v5)
        }
        else {
          next = right[current];
          weight1[length1] = wr[current];
          weight2[length2] = wr[current];
        }
        tmp = current;
        current = next;

      } // End while

// Add vertex i at the end of path1 and path2. Also add edge (i,v1) at the end of path2
      length1++;
      length2++;
      path1[length1] = current;
      path2[length2] = current;

      if (right[current] == tmp) { // Moving next forward to v5
        weight2[length2] = wl[current];  // Picking up weight of edge (v4,v5)
        next = left[current];
      }
      else {
        weight2[length2] = wr[current];
        next = right[current];
      }
      length2++;
      path2[length2] = next;

// Now pick the best of path1 and path2
      if (dyn_prog(length1,path1,weight1,match) > dyn_prog(length2,path2,weight2,match))
        dyn_prog(length1,path1,weight1,match);
      else
        dyn_prog(length2,path2,weight2,match);

      continue;
    } // End treatment of cycle

    if ((right[i] != 0) && (left[i] !=0)) // This is not an endpoint, skip it
      continue;

    if ((right[i] == 0) && (left[i] ==0)) // This is a singleton, skip it
      continue;

    if (other[i] == -2) // This is a path that has been handled from the other endpoint
      continue;

// We now have a path starting in i, must pick up every vertex and mark other endpoint with -2
    current = i;
    if (right[i] == 0) {
      next = left[i];
      weight1[0] = wl[i];
    }
    else {
      next = right[i];
      weight1[0] = wr[i];
    }
    length1 = 0;
    while (current != other[i]) { // Move through the path picking up vertices and edge weights
      path1[length1] = current;
      length1++;
      int tmp = next;
      if (right[next] != current) {
        next = right[next];
        weight1[length1] = wr[tmp];
      } 
      else {
        next = left[next];
        weight1[length1] = wl[tmp];
      }
      current = tmp;
    }
    path1[length1] = current; // Add the final vertex
    other[current] = -2;     // Mark last vertex so we only traverse path in one direction

    dyn_prog(length1,path1,weight1,match);

  }

// Calculate the cost of the matching

  float glob_sum = 0.0;
  int k;
  for(i=1;i<=n;i++) {   
    if ((match[i] != 0) && (i < match[i]))
      for(k=ver[i];k<ver[i+1];k++)
        if (edges[k] == match[i]) {
          glob_sum += weight[k];
          if (match[edges[k]] != i) {
            printf("**** Error in GPA! \n");
            printf("%d is matched to %d but %d is matched to %d \n",i,match[i],match[i],match[edges[k]]);
            return;
          }
        }
  }

//  printf("GPA: Dynamic programming gave matching of weight %f \n",glob_sum);

  for(i=1;i<=n;i++) {
    if (match[i] == 0) {
      float best = -1.0;
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

//  printf("GPA: After removing single nodes weight is %f \n",glob_sum);

}
