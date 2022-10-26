// Row sums of matrix
vector row_sums(matrix X) {
  vector[rows(X)] s = X * rep_vector(1, cols(X));
  return s;
}
