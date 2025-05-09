// R-style conditional subsetting
array[] int rsub(array[] int x, array[] int cond) {
  array[sum(cond)] int xsub;
  int pos;
  pos = 1;
  for (i in 1:size(x))
    if (cond[i])
    {
      xsub[pos] = x[i];
      pos = pos + 1;
    }
  return(xsub);
}
