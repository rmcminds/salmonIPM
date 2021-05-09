// Quantiles of a vector
real quantile(vector v, real p) {
  int N = num_elements(v);
  real Np = round(N*p);
  real q;
  
  for(i in 1:N) {
    if(i - Np == 0.0) q = v[i];
  }
  return(q);
}
