// Vectorized logical &&
array[] int vand(array[] int cond1, array[] int cond2) {
  array[size(cond1)] int cond1_and_cond2;
  for(i in 1:size(cond1)) cond1_and_cond2[i] = cond1[i] && cond2[i];
  return(cond1_and_cond2);
}
