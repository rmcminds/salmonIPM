// Vectorized logical equality
array[] int veq(array[] int x, int y) {
  int xeqy[size(x)];
  for(i in 1:size(x))
    xeqy[i] = x[i] == y;
  return(xeqy);
}
