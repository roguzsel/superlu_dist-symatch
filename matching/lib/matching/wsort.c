// 
// Sorting neighbor lists by weight using insertion sort
//

void wsort(int n,int *ver,int *edges,double *weight,int **s_edges,double **s_weight) {

  int i,j;
  int cur,e;
  double we;
  void quick();

// Allocate space for sorted edge list

  (*s_edges) = (int *) malloc(sizeof(int)*(ver[n+1]));
  (*s_weight) = (double *) malloc(sizeof(double)*(ver[n+1]));

  for(i=0;i<ver[n+1];i++) {                // Copy edges and weights into s_edges and s_weight
    (*s_edges)[i] = edges[i]; 
    (*s_weight)[i] = weight[i]; 
  }
 
  for(i=1;i<=n;i++) {                // Loop over vertices
    for(j=ver[i]+1;j<ver[i+1];j++) { // Loop over neighbors of vertex i, except for first one

/*
    quick(edges,weight,ver[i],ver[i+1]-1);
*/
      cur = j;
      we = (*s_weight)[cur];
      e = (*s_edges)[cur];
      while ((we > (*s_weight)[cur-1]) || ((we == (*s_weight)[cur-1]) && (e < (*s_edges)[cur-1]))) {  // Swap cur with cur-1
        (*s_weight)[cur] = (*s_weight)[cur-1];
        (*s_edges)[cur] = (*s_edges)[cur-1];
        cur--;
        if (cur == ver[i])
          break;
      } // loop until neighbor is in the right place
      (*s_weight)[cur] = we;
      (*s_edges)[cur] = e;
    } // loop over neighbors
  } // loop over vertices

}

void swap(int array[],double weight[], int x, int y) {
  
  int a = array[x];
  array[x] = array[y];
  array[y] = a;

  double b = weight[x];
  weight[x] = weight[y];
  weight[y] = b;
}

void quick(int array[], double weight[], int start, int end){


    if(start < end){
        int l=start+1, r=end; 
        double p = weight[start];
        while(l<r){
            if(weight[l] <= p)
                l++;
            else if(weight[r] >= p)
                r--;
            else
                swap(array,weight,l,r);
        }
        if(weight[l] < p){
            swap(array,weight,l,start);
            l--;
        }
        else{
            l--;
            swap(array,weight,l,start);
        }
        quick(array, weight, start, l);
        quick(array, weight, r, end);
    }
}

