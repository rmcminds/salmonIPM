// R-style conditional subsetting
int[] rsub(int[] x, int[] cond) {
  int xsub[sum(cond)];
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
