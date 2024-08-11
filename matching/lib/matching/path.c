// Path growing algorithm with DP
//

void path(int n,int *ver,int *edges,double *weight,int *match ) {

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

// Start of matching algorithm

//  printf("Starting sequential matching algorithm \n");

  pw = (double *) malloc((n+1)* sizeof(double));
  if (pw == NULL) {
    printf("Unable to allocate space for pw in path() \n");
    return;
  }
  path = (int *) malloc((n+1)* sizeof(int));
  if (path == NULL) {
    printf("Unable to allocate space for path in path() \n");
    return;
  }
  visited = (int *) malloc((n+1)* sizeof(int));
  if (visited == NULL) {
    printf("Unable to allocate space for visited in path() \n");
    return;
  }

  for(i=1;i<=n;i++) {   // Set that every node is unmatched
    visited[i] = false;
    match[i] = 0;
  }


  for(i=1;i<=n;i++) {                 // Loop over vertices


    if (visited[i])
      continue;
    int done = false;
    int next = i;		// Start a new path
    length = -1;
    while (!done) {
      length++;

      path[length] = next;   	// Store the edge of the path
      visited[next] = true; 
      heaviest = -1.0;
      for(j=ver[next];j<ver[next+1];j++) { // Loop over neighbors of vertex next, find heaviest available edge
        y = edges[j];
        if ((!visited[y]) && (weight[j]>heaviest)) {     // Check if vertex is unmatched and if it gives a new heaviest edge
          heaviest = weight[j];
          best = y;
        }
      }

      if (heaviest < 0.0)  	// If there was no available neighbor then terminate
        done = true;
      else {			// Store heaviest edge and move on
        next = best;
        pw[length] = heaviest;
      }
    } // End while

// Now do dynamic programming
    sdyn_prog(length,path,pw,match);

  } // loop over vertices


  int k;

/*
  for(i=1;i<=n;i++) {   
    if ((match[i] != 0) && (i < match[i]))
      for(k=ver[i];k<ver[i+1];k++)
        if (edges[k] == match[i]) {
          glob_sum += weight[k];
          if (match[edges[k]] != i) {
            printf("*** Error in path growing \n");
            printf("%d is matched with %d, but %d is matched with %d \n",i,match[i],match[i],match[match[i]]);
          }
        }
  }
*/

// Now do greedy matching on any remaining free vertices

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

//  printf("Have matched %d vertices out of %d \n",nr_matched,n);
}


double sdyn_prog(int l1,int *path1,double *weight1,int *match) {

  double temp;
  double w_old;
  double w_new;
  int i,k;
//  int in[l1];
  int *in;
  int partner;

  w_old = 0.0;
  w_new = weight1[0];

  in = (int *) malloc((l1)* sizeof(int)); 
  if (in == NULL) {
    printf("Unable to allocate space for in in dyn_prog \n");
    return(0.0);
  } 
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
  free(in);
  return w_new;  // Return cost of best matching
}


