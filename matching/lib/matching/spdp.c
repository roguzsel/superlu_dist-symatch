// weighted parallel algorithm
//
// Two round suitor based matching followed by dynamic programming
// Assuming sorted neighbor lists

void spdp(int n,int *ver,int *edges,int *s,double *ws,int *s2, double *ws2,double *weight, int *used,int *match,omp_lock_t *nlocks,int *par,int *path1,double *weight1) {


  double cost_matching();

  int i,j,e;
  int partner;     // Prospective candidate to match with   
  int x,y, done, next_vertex;
  int count,current;
  // int *path1;
  int l1,l2;
  // double *weight1;
  double cum_weight = 0.0;
  double dyn_prog();

  int my_id = omp_get_thread_num();

// Allocate these outside
//  path1 = (int *) malloc(n * sizeof(int));
//  weight1 = (double *) malloc(n * sizeof(double));



// Start of matching algorithm
#pragma omp for schedule(static) private(i)
  for(i=0;i<=n;i++) {   
    s[i] = 0;           // Set that no node is trying to match with i 
    ws[i] = 0.0;        // Set the current weight of best suitor edge to 0

    s2[i] = 0;           // Set that no node is trying to match with i 
    ws2[i] = 0.0;        // Set the current weight of best suitor edge to 0
    omp_init_lock(&nlocks[i]); // Initialize locks

    used[i] = false;     // Set that node i has not been used in the dynamic programming part
    match[i] = 0;
    par[i] = ver[i];     // Set where in edge list to start searching for a partner 
  }
//  count = 0;


  int *next;
  next = par;

#pragma omp barrier

// First round of matching, this is regular suitor based parallel matching

// #pragma omp for schedule(dynamic,700) private(x)
#pragma omp for schedule(static) private(x)
  for(x=1;x<=n;x++) {                // Loop over vertices
    done = false;
    i = x;
    while (!done) {
/*
      double heaviest = ws[i];   // Start search from current suitor
      partner = s[i];
*/
      y = next[i];          // y points to the next possible neighbor (in the edge list of i) that i can match with
      e = edges[y];         // Get the id of the neighbor that we will try to match with

// Stop if there are no more edges or if we found a candidate
/*
      while ((y < ver[i+1]) && ((weight[y] < ws[e]) || ((weight[y] == ws[e]) && (i > s[e])))) {
        y++;              // Move to next position in edge list
        e = edges[y];     // Get id of neighbor
      }
*/
      while ((y < ver[i+1]) && (weight[y] < ws[e])) {
        y++;              // Move to next position in edge list
        e = edges[y];     // Get id of neighbor
      }


      done = true;
      if (y < ver[i+1]) {            // Test if we found a possible partner
        next[i] = y+1;               // Set where to search from next time i needs a partner
        omp_set_lock(&nlocks[e]);    // Locking e 

        if ((weight[y] > ws[e]) || ((weight[y] == ws[e]) && (i > s[e]))) {
          if (s[e] != 0) {             // True if e already had a suitor
            next_vertex = s[e];        // Pick up the old suitor
            done = false;
          }
          s[e] = i;                  // i is now the current suitor of e
          ws[e] = weight[y];           // Store the weight of (i,e) as the weight given by the current suitor for e
        }
        else {   // can no longer use the result for node i, must find new partner
          done = false;
          next_vertex = i;
        }
        omp_unset_lock(&nlocks[e]);  // Unlocking e 
      } // if y
      if (!done) {  // Continue with the next vertex
        i = next_vertex;
      }
    } // while not done
  } // loop over vertices


#pragma omp barrier

// Second round of matching



// #pragma omp for schedule(dynamic,700) private(x)
#pragma omp for schedule(static) private(x)
  for(x=1;x<=n;x++) {                // Loop over vertices
    next[x] = ver[x];     // Set where in edge list to start searching for a partner 
    done = false;
    i = x;
    while (!done) {
/*
      double heaviest = ws[i];   // Start search from current suitor
      partner = s[i];
*/
      y = next[i];          // y points to the next possible neighbor (in the edge list of i) that i can match with
      e = edges[y];         // Get the id of the neighbor that we will try to match with

// Stop if there are no more edges or if we found a candidate
/*
      while ((y < ver[i+1]) && ((weight[y] < ws2[e]) || ((weight[y] == ws2[e]) && (i > s2[e])) || (y == s[i]))) {
        y++;              // Move to next position in edge list
        e = edges[y];     // Get id of neighbor
      }
*/
      while ((y < ver[i+1]) && ((weight[y] < ws2[e]) || (e == s[i]))) {
        y++;              // Move to next position in edge list
        e = edges[y];     // Get id of neighbor
      }


      done = true;
      if (y < ver[i+1]) {            // Test if we found a possible partner
        next[i] = y+1;               // Set where to search from next time i needs a partner
        omp_set_lock(&nlocks[e]);    // Locking e 

        if ((weight[y] > ws2[e]) || ((weight[y] == ws2[e]) && (i > s2[e]))) {
          if (s2[e] != 0) {             // True if e already had a suitor
            next_vertex = s2[e];        // Pick up the old suitor
            done = false;
          }
          s2[e] = i;                  // i is now the current suitor of e
          ws2[e] = weight[y];           // Store the weight of (i,e) as the weight given by the current suitor for e
        }
        else {   // can no longer use the result for node i, must find new partner
          done = false;
          next_vertex = i;
        }
        omp_unset_lock(&nlocks[e]);  // Unlocking e 
      } // if y
      if (!done) {  // Continue with the next vertex
        i = next_vertex;
      }
    } // while not done
  } // loop over vertices



#pragma omp barrier
/*
#pragma omp master 
 { double cost = cost_matching(n,ver,edges,weight,s2);
   printf("Cost of second matching is %4.2f \n",cost);
 }
*/

// Starting dynamic programming



// #pragma omp for schedule(static) private(i)
// #pragma omp for schedule(dynamic,700) 
#pragma omp for schedule(static) 
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path

    par[i] = i;
    if ((s[i] != 0) && (s2[i] != 0)) continue;  // If this vertex is not an endpoint then move on

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

      if (path1[0] < path1[l1]) {  // Make sure DP only starts from one end of the path
        // printf("%d, Doing DP on path of length %d \n",my_id,l1);
      
        dyn_prog(l1,path1,weight1,match);
      }

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

      if (path1[0] < path1[l1]) {  // Only start DP from one end of the path
        // printf("%d: Doing DP on path of length %d \n",my_id,l1);
     
        dyn_prog(l1,path1,weight1,match);
      }

    }
    else if ((s2[i] == 0) && (s[i] == 0)) {  // Found single vertex
      used[i] = true;
    }
  }


#pragma omp barrier

// Now do the cycles


  int *left;
  left = (int *) malloc((1+n/omp_get_num_threads()) * sizeof(int));  // Space for holding one endpoint of the level 2 edges

  if (left == NULL) {
    printf("Unable to allocate space for left in pdp \n");
    return;
  }

  int remaining = 0;

// Do union on every level one edge and store one end point of the level 2 edges

// #pragma omp for schedule(dynamic,700) private(i)
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (used[i]) continue;  // If this vertex has been used previously then skip it

    if (i > s[i]) {   // Link i and s[i]
      par[i] = i;
      par[s[i]] = i;
    }

    if (i < s2[i]) {  // Store i
      left[remaining] = i;
      remaining++;
    }
  }


#pragma omp barrier

// Run Rem's algorithm on every level two cycle edge, this will merge the cycles into sets

  for(i=0;i<remaining;i++) {  
    x = left[i];
    y = s2[x];

    while (par[x] != par[y]) {   // Check if x and y have the same parent
      if (par[x] < par[y]) {     // Find the one with the smaller parent
        if (x == par[x]) {       // If x is a root we can link

          omp_set_lock(&nlocks[x]); // locking for x
          int p_set = false;   // Check if someone else got to x first
          if (x == par[x]) {   // if x is still a root
            par[x] = par[y];   // linking
            p_set = true;      // indicate success
//              nr_union++;
          }
          omp_unset_lock(&nlocks[x]); // un-locking for x

          if (p_set)
            break;
        }

        int z = par[x];         // Splicing
        par[x] = par[y];
        x = z;
      }
      else {
        if (y == par[y]) {      // If y is a root we can link

          omp_set_lock(&nlocks[y]); // locking for y
          int p_set = false;   // Check if someone else got to x first
          if (y == par[y]) {   // if y is still a root
            par[y] = par[x];   // linking
            p_set = true;      // indicate success
//            nr_union++;
          }
          omp_unset_lock(&nlocks[y]); // un-locking for y

          if (p_set)
            break;
        } // if

        int z = par[y];         // Splicing
        par[y] = par[x];
        y = z;

      } // else
    }  // while loop
  } // loop over edges


#pragma omp barrier

// Now look for the cycles, each one is represented by its root vertex

  int nr_cycle = 0;
// #pragma omp for schedule(dynamic,700) private(i)
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (i != par[i]) continue;  // Only use root vertices as starting points of a cycle

    nr_cycle++;
    current = i;
    l1 = 0;

// Ready to start loop
      
    while (!used[current]) {   // Process two edges (vi,v(i+1)) and (v(i+1),v(i+2)), at a time until we reach the start of the cycle

// Put vertex v_i and edge (vi,v(i+1)) in the path

      used[current] = true;    		// Set vi as used
      path1[l1] = current;		// Put vi in path1
      weight1[l1] = ws[current];	// Put weight of (vi,v(i+1)) in the path

      current = s[current];    		// Move to v(i+1)

      l1++;				// Increase number of vertices
      used[current] = true;    		// Set v(i+1) as used
      path1[l1] = current;		// Put v(i+1) in path1
      weight1[l1] = ws2[current];	// Put weight of (v(i+1),v(i+2) in the path

      current = s2[current]; 		// Move to v(i+2)
      l1++;
    } // End while loop

    path1[l1] = path1[0];		// Put the first vertex last to complete the cycle


// Now do dynamic programming on the cycle 
    cycle_dyn_prog(l1,path1,weight1,match);
 
  } // End of loop over vertices that looks for unprocessed cycles


#pragma omp barrier

// Now pick up the remaining singletons and match these greedily

  int k;
  remaining = 0;
// #pragma omp for schedule(dynamic,700) private(i)
#pragma omp for schedule(static) private(i)
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
        left[remaining] = i;    // Store the vertices that need to be checked
        remaining++;
        left[remaining] = partner;    // Store the vertices that need to be checked
        remaining++;
      }

    }
  }


#pragma omp barrier

// Loop through the newly matched vertices and unmatch any vertex where there was a conflict

  for(i=0;i<remaining;i++) {
    if (match[match[left[i]]] != left[i]) {
      match[left[i]] = 0;
    }
  }

  free(left);
//  free(path1);
//  free(weight1);

#pragma omp barrier
  return;  


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

