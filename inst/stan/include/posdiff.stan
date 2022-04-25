// Constrain a difference to be positive by adjusting the subtrahend if necessary
// Similar to ADMB posfun() but returns adjusted subtrahend instead of difference 
// For posfun() see https://github.com/kaskr/adcomp/issues/7
real posdiff(real x, real y, real eps) {
  real d = x - y;
  if(d > 0) return(y);
  else return(x - eps/(2 - d/eps));
}
