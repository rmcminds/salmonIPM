// Vectorized logical equality
array[] int veq(array[] int x, int y) {
  array[size(x)] int xeqy;
  for(i in 1:size(x)) xeqy[i] = x[i] == y;
  return(xeqy);
}
