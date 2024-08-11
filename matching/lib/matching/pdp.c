// weighted parallel algorithm
//
// Two round suitor based matching followed by dynamic programming

void pdp(int n,int *ver,int *edges,int *s,double *ws,int *s2, double *ws2,double *weight, int *used,int *match,omp_lock_t *nlocks,int *par,int *path1,double *weight1) {


  double cost_matching();

  int i,j;
  int partner;     // Prospective candidate to match with   
  int x,y, done, next_vertex;
  int count,current;
  // int *path1;
  int l1,l2;
  // double *weight1;
  double cum_weight = 0.0;
  double dyn_prog();

  int my_id = omp_get_thread_num();

  double mt1,mt2;


  // path1 = (int *) malloc(n * sizeof(int));
  // weight1 = (double *) malloc(n * sizeof(double));

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
    par[i] = 0;          // Parrent pointer in Rem's algorithm
  }


#pragma omp barrier

// First round of matching, this is regular suitor based parallel matching
/*
  #pragma omp master
        { mt1 = omp_get_wtime();
        }
*/

//#pragma omp for schedule(dynamic,700) private(x)
 #pragma omp for schedule(static) private(x)
  for(x=1;x<=n;x++) {                // Loop over vertices
    i = x;
    done = false;
    while (!done) {

/*
      omp_set_lock(&nlocks[i]);    // Locking i 
      double heaviest = ws[i];   // Start search from current suitor
      partner = s[i];
      omp_unset_lock(&nlocks[i]); // Unlocking i 
*/
      double heaviest = -1.0;

      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int y = edges[j];            // y is the current neighbor of i
/*
        if ((weight[j]>heaviest) && (ws[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y
*/


        if ((weight[j] < heaviest) || (weight[j] < ws[y])) // If the weight of the current edge is lighter than the current best or
          continue;                                        // if we cannot offer y something better, move on.

        if ((weight[j] == heaviest) && (y > partner))      // When equal weight, give priority to highest index
          continue;

        if ((weight[j] == ws[y]) && (i > s[y]))            // Must have higher index to beat previous suitor when equal weight
          continue;

// Check if w(i,y) is the best so far, and if it is a better option for y
/*
        if (((weight[j] > heaviest) || ((weight[j] == heaviest) && (y > partner))) &&
            ((weight[j] > ws[y]) || ((weight[j] == ws[y]) && (i > s[y])))) {
*/
        heaviest = weight[j];      // Store the best so far
        partner = y;

      } // loop over neighbors

      done = true;
      if (heaviest > 0) {
        omp_set_lock(&nlocks[partner]);    // Locking partner

        // if ((heaviest > ws[partner]) || ((heaviest == ws[partner]) && (i>s[partner]))) {
        if (heaviest >= ws[partner]) { // Must have >= because i might start from its best suitor
          if (s[partner] != 0) {
            next_vertex = s[partner];
            done = false;
          }
          s[partner] = i;          // i is now the current suitor of s
          ws[partner] = heaviest;  // The weight of the edge (i,partner) 
        }
        else {   // can no longer use the result for node i, must find new partner
          done = false;
          next_vertex = i;
        }
        omp_unset_lock(&nlocks[partner]); // Unlocking partner
      }
      if (!done) {  // Continue with the next vertex
        i = next_vertex;
      }
    }  // while not done
  } // loop over vertices

#pragma omp barrier

/*
  #pragma omp master
        { mt2 = omp_get_wtime();
          printf("First round took %lf \n",mt2-mt1);
        }
*/
// Second round of matching

// #pragma omp for schedule(dynamic,700) private(x)
#pragma omp for schedule(static) private(x)
  for(x=1;x<=n;x++) {                // Loop over vertices
    i = x;
    done = false;
    while (!done) {

/*
      omp_set_lock(&nlocks[i]);    // Locking i 
      double heaviest = ws2[i];   // Start search from current suitor
      partner = s2[i];
      omp_unset_lock(&nlocks[i]); // Unlocking i 
*/


      double heaviest = -1.0;
      partner = 0;
      for(j=ver[i];j<ver[i+1];j++) { // Loop over neighbors of vertex i
        int y = edges[j];            // y is the current neighbor of i

/*
        if ((weight[j]>heaviest) && (ws2[y]<weight[j]) && (y != s[current])) {
          heaviest = weight[j];      
          partner = y;              
      } 
*/
//        if ((weight[j]>heaviest) && (ws2[y]<weight[j])) {// Check if w(i,y) is the best so far, and if it is a better option for y

        if ((weight[j] < heaviest) || (weight[j] < ws2[y])) // If the weight of the current edge is lighter than the current best or
          continue;                                         // if we cannot offer y something better, move on.

        if (y == s[i])                                      // If this edge was used in the first matching then skip it
          continue;                                     

        if ((weight[j] == ws2[y]) && (i > s2[y]))           // Must have lower index to beat previous suitor when equal weight
          continue;

        if ((weight[j] == heaviest) && (y > partner))       // When equal weight, give priority to lowest index
          continue;



// Check if w(i,y) is the best so far, and if it is a better option for y

//        if (((weight[j] > heaviest) || ((weight[j] == heaviest) && (y > partner))) &&
//            ((weight[j] > ws2[y]) || ((weight[j] == ws2[y]) && (i > s2[y])))) {

        heaviest = weight[j];      // Store the best so far
        partner = y;

      } // loop over neighbors

      done = true;
      if (heaviest > 0) {
        omp_set_lock(&nlocks[partner]);    // Locking partner

        // if ((heaviest > ws2[partner]) || ((heaviest == ws2[partner]) && (i>s2[partner]))) {
        if (heaviest >= ws2[partner]) { // Must have >= because i might start from its best suitor
          if (s2[partner] != 0) {
            next_vertex = s2[partner];
            done = false;
          }
          s2[partner] = i;          // i is now the current suitor of s
          ws2[partner] = heaviest;  // The weight of the edge (i,partner) 
        }
        else {   // can no longer use the result for node i, must find new partner
          done = false;
          next_vertex = i;
        }
        omp_unset_lock(&nlocks[partner]); // Unlocking partner
      }
      if (!done) {  // Continue with the next vertex
        i = next_vertex;
      }
    }  // while not done
  } // loop over vertices


#pragma omp barrier

// Starting dynamic programming
/*
  #pragma omp master
        { mt1 = omp_get_wtime();
          printf("Second round took %lf \n",mt1-mt2);
        }
*/

// #pragma omp for schedule(dynamic,700) private(i)
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
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
        dyn_prog(l1,path1,weight1,match);
      }

    }
    else if ((s2[i] == 0) && (s[i] == 0)) {  // Found single vertex
      used[i] = true;
    }
  }

/*
  #pragma omp master
        { mt2 = omp_get_wtime();
          printf("Dynamic programming took %lf \n",mt2-mt1);
        }
*/

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

// #pragma omp for schedule(dynamic,700) private(i)
#pragma omp for schedule(static) private(i)
  for(i=1;i<=n;i++) {  // Check each vertex as a starting point of a path
    if (i != par[i]) continue;  // Only use root vertices as starting points of a cycle

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
//#pragma omp for schedule(dynamic,700) private(i)
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
  // free(path1);
  // free(weight1);

  return;  


}

