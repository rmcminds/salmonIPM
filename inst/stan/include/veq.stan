// Vectorized logical equality
int[] veq(int[] x, int y) {
  int xeqy[size(x)];
  for(i in 1:size(x))
    xeqy[i] = x[i] == y;
  return(xeqy);
}
