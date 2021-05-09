// Left multiply vector by matrix
// works even if size is zero
vector mat_lmult(matrix X, vector v)
{
  vector[rows(X)] Xv;
  Xv = rows_dot_product(X, rep_matrix(to_row_vector(v), rows(X)));
  return(Xv);
}
