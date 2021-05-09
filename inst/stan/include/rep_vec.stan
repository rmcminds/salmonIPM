// Return the size m*n column vector consisting of m copies of x,
// where x is a column vector of size n
vector rep_vec(vector x, int m) {
  vector[num_elements(x)*m] y;
  y = to_vector(rep_matrix(x,m));
  return(y);
}
