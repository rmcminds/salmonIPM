// spawner-recruit functions
real SR(int SR_fun, real alpha, real Rmax, real S, real A) {
  real R;
  
  if(SR_fun == 1)      // discrete exponential
    R = alpha*S;
  else if(SR_fun == 2) // Beverton-Holt
    R = alpha*S/(1 + alpha*S/(A*Rmax));
  else if(SR_fun == 3) // Ricker
    R = alpha*S*exp(-alpha*S/(A*e()*Rmax));
  
  return(R);
}
