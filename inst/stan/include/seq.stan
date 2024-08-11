// Equivalent of R operator ":"
array[] int seq(int from, int to) {
  int x[to-from+1];
  for(i in from:to) x[i-from+1] = i;
  return(x);
}
