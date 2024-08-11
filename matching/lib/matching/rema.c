// 2-augmenting code
//
// This code assumes that s contains a matching with ws giving the weight
//

void rema(int n,int *ver,int *edges,int *s,double *ws,double *weight,int *p1) {

  return;
}


void addEdge(int *s,double *ws,int x,int y,double newW,int *wQ,int *round,int *wLast,int curRound) {
   // printf("Replacing (%d,%d) and (%d,%d)=(%d,%d) \n",x,s[x],y,s[y],s[y],s[s[y]]);
   // printf("of weight %lf=%lf  and %lf=%lf \n",ws[x],ws[s[x]],ws[y],ws[s[y]]);
   // printf("Matching %d and %d of weight %lf \n",x,y,newW);
  s[s[x]] = 0;
  s[s[y]] = 0;
  ws[s[x]] = 0.0;
  ws[s[y]] = 0.0;
  int a = s[x];
  int b = s[y];
  s[x] = y;
  s[y] = x;
  ws[x] = newW;
  ws[y] = newW;
/*
  if ((a != 0) && (round[a] <= curRound)) {
    round[a] = curRound+1;
    wQ[*wLast] = a;
    *wLast = *wLast + 1;
  }
  if ((b != 0) && (round[b] <= curRound)) {
    round[b] = curRound+1;
    wQ[*wLast] = b;
    *wLast = *wLast + 1;
  }
*/
  // (Load up the edge (x,y) to be checked in the next round, unless it has already been put in the queue
  if ((round[x]>curRound) || (round[y]>curRound))
    return;
  round[x] = curRound+1;
  wQ[*wLast] = x;           // Put x in the queue
  *wLast = *wLast + 1;
//  printf("Added %d to queue \n",x);
}
