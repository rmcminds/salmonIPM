// Generalized normal (aka power-exponential) unnormalized log-probability
real pexp_lpdf(real y, real mu, real sigma_R, real shape) {
  return(-(fabs(y - mu)/sigma_R)^shape);
}
