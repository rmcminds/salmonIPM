// Column sums of matrix
row_vector col_sums(matrix X) {
  row_vector[cols(X)] s = rep_row_vector(1, rows(X)) * X;
  return s;
}
