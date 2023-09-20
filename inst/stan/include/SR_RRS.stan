// spawner-recruit functions
real SR(int SR_fun, int[] RRS, real alpha, real alpha_W, real alpha_H, 
        real Rmax, real Rmax_W, real Rmax_H, real S, real S_W, real S_H, real A) {
  real R;
  
  // no H vs. W differences
  if(!RRS[1] && !RRS[2])
  {
    if(SR_fun == 1) // discrete exponential
      R = alpha*S;
    if(SR_fun == 2) // Beverton-Holt
      R = alpha*S/(1 + alpha*S/(A*Rmax));
    if(SR_fun == 3) // Ricker
      R = alpha*S*exp(-alpha*S/(e()*A*Rmax));
  }

  // intrinsic productivity: H != W, max production: H == W
  if(RRS[1] && !RRS[2])
  {
    if(SR_fun == 1) // discrete exponential
      R = alpha_W*S_W + alpha_H*S_H;
    if(SR_fun == 2) // Beverton-Holt (Leslie-Gower)
      R = (alpha_W*S_W + alpha_H*S_H) / (1 + (alpha_W*S_W + alpha_H*S_H)/(A*Rmax));
    if(SR_fun == 3) { // Ricker
      real log_alpha_W = log(alpha_W);
      real log_alpha_H = log(alpha_H);
      real DD = 1 - alpha_W*S_W/(e()*A*Rmax*log_alpha_W) - alpha_H*S_H/(e()*A*Rmax*log_alpha_H);
      R = S_W*exp(log_alpha_W*DD) + S_H*exp(log_alpha_H*DD);
    }
  }

  // intrinsic productivity: H == W, max production: H != W
  if(!RRS[1] && RRS[2])
  {
    if(SR_fun == 1) // discrete exponential
      R = alpha*S;
    if(SR_fun == 2) // Beverton-Holt (Leslie-Gower)
      R = alpha*S / (1 + alpha*S_W/(A*Rmax_W) + alpha*S_H/(A*Rmax_H));
    if(SR_fun == 3) // Ricker
      R = alpha*S * exp(-alpha*S_W/(e()*A*Rmax_W) - alpha*S_H/(e()*A*Rmax_H));
  }

  // intrinsic productivity and max production: H != W
  if(RRS[1] && RRS[2])
  {
    if(SR_fun == 1) // discrete exponential
      R = alpha_W*S_W + alpha_H*S_H;
    if(SR_fun == 2) // Beverton-Holt (Leslie-Gower)
      R = (alpha_W*S_W + alpha_H*S_H) / (1 + alpha_W*S_W/(A*Rmax_W) + alpha_H*S_H/(A*Rmax_H));
    if(SR_fun == 3) { // Ricker
      real log_alpha_W = log(alpha_W);
      real log_alpha_H = log(alpha_H);
      real DD = 1 - alpha_W*S_W/(e()*A*Rmax_W*log_alpha_W) - alpha_H*S_H/(e()*A*Rmax_H*log_alpha_H);
      R = S_W*exp(log_alpha_W*DD) + S_H*exp(log_alpha_H*DD);
    }
  }
  
  return(R);
}
