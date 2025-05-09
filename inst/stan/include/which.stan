// Equivalent of R: which(cond), where sum(cond) == 1
int which(array[] int cond) {
  int which_cond;
  for(i in 1:size(cond)) if(cond[i]) which_cond = i;
  return(which_cond);
}
