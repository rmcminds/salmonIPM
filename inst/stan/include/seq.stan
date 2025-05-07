// Equivalent of R operator ":"
array[] int seq(int from, int to) {
  array[to-from+1] int x;
  for(i in from:to) x[i-from+1] = i;
  return(x);
}
