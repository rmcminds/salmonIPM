functions {
  #include /include/SR.stan
  #include /include/gnormal_lpdf_vec.stan
}

data {
  int<lower=1> SR_fun;          // S-R model: 1 = exponential, 2 = BH, 3 = Ricker, 4 = Hassell
  int<lower=1> N;               // total number of cases in all pops and years
  array[N] int<lower=1,upper=N> pop;  // population index
  array[N] int<lower=1,upper=N> year; // brood year index
  int<lower=1,upper=N> N_fit;   // number of cases used in fitting (non-missing S_obs and R_obs)
  array[N_fit] int<lower=1,upper=N> which_fit; // cases used in fitting
  vector<lower=0>[N] S_obs;     // observed annual total spawner abundance (not density)
  vector<lower=0>[N] R_obs;     // total natural recruit abundance (not density), including harvest and broodstock removals
  vector[N] A;                  // habitat area associated with each spawner abundance obs
  array[N] int<lower=0,upper=1> S_NA; // logical indicating whether S_obs is missing and should be simulated
  array[N] int<lower=0,upper=1> R_NA; // logical indicating whether R_obs is missing and should be simulated
  int<lower=2> N_age;           // number of adult age classes
  int<lower=2> max_age;         // maximum adult age
  matrix<lower=0,upper=1>[max(pop),N_age] p_pop_obs;  // average recruit age distributions for each pop 
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  array[N_age] int<lower=2> ages;          // adult ages
  real mu_Rmax = max(log(R_obs[which_fit] ./ A[which_fit]));  // prior log-mean of Rmax
  real sigma_Rmax = sd(log(R_obs[which_fit] ./ A[which_fit]));  // prior log-SD of Rmax
  
  for(a in 1:N_age)
    ages[a] = max_age - N_age + a;
}

parameters {
  vector<lower=0>[N_pop] alpha;           // intrinsic productivity
  vector<lower=0>[N_pop] Rmax;            // maximum recruitment
  vector<lower=-1,upper=1>[N_pop] rho_R;  // AR(1) coefs of recruitmentresiduals
  vector<lower=0>[N_pop] sigma_R;         // recruitment residual error SD
  vector<lower=0>[SR_fun == 4 ? 1 : 0] b; // Hassell shape parameter
}

transformed parameters {
  vector<lower=0>[N] R_hat;     // expected recruit abundance (not density) by brood year
  vector<lower=0>[N] R_ar1;     // expected recruit abundance, taking AR(1) errors into account
  vector<lower=0>[N] sigma_ar1; // residual error SD for each observation
  
  // Predict recruitment
  R_hat = rep_vector(0,N);
  R_ar1 = rep_vector(0,N);
  sigma_ar1 = rep_vector(0,N);
    
  for(i in 1:N_fit)
  {
    int ii = which_fit[i];
    int il = which_fit[i-1];
    
    R_hat[ii] = SR(SR_fun, {0,0}, alpha[pop[ii]], 0, 0, Rmax[pop[ii]], 0, 0, 
                   S_obs[ii], 0, 0, A[ii], b);
    if(i==1 || pop[il] != pop[ii])
    {
      R_ar1[ii] = R_hat[ii];
      sigma_ar1[ii] = sigma_R[pop[ii]]/sqrt(1 - rho_R[pop[ii]]^2);
    }
    else
    {
      real err;    // residual at the last non-missing observation
      int dt;      // number of years since last non-missing observation
      real rho2j;  // sum of powers of rho_R
      
      err = log(R_obs[il]) - log(R_hat[il]);
      dt = ii - il;
      R_ar1[ii] = R_hat[ii]*exp(rho_R[pop[ii]]^dt * err);
      
      rho2j = 0;
      for(j in 0:(dt-1))
        rho2j = rho2j + rho_R[pop[ii]] ^ (2*j);
      sigma_ar1[ii] = sigma_R[pop[ii]] * sqrt(rho2j);
    }
  }
}

model {
  // Priors
  alpha ~ lognormal(2.0,2.0);
  Rmax ~ lognormal(mu_Rmax, sigma_Rmax);
  rho_R ~ gnormal(0,0.85,20);   // mildly regularize rho_R to ensure stationarity
  sigma_R ~ normal(0,3);
  b ~ exponential(1);

  // Likelihood
  R_obs[which_fit] ~ lognormal(log(R_ar1[which_fit]), sigma_ar1[which_fit]);
}

generated quantities {
  vector[N] S_sim;    // simulated spawners
  vector[N] R_sim;    // simulated recruits
  vector[N] err_sim;  // simulated AR(1) residual errors
  
  S_sim = S_obs;
  R_sim = R_obs;
  err_sim = rep_vector(0,N);
  
  for(i in 1:N)
  {
    if(S_NA[i] == 1)
    {
      if(i >= max_age && pop[i-max_age] == pop[i])
      {
        S_sim[i] = 0;
        for(a in 1:N_age)
          S_sim[i] = S_sim[i] + R_sim[i-ages[a]]*p_pop_obs[pop[i],a];
      }
    }
    
    if(i == 1 || pop[i-1] != pop[i])
      err_sim[i] = normal_rng(0, sigma_R[pop[i]]);
    else
      err_sim[i] = normal_rng(rho_R[pop[i]]*err_sim[i-1], sigma_R[pop[i]]);
    
    if(R_NA[i] == 1)
      R_sim[i] = SR(SR_fun, {0,0}, alpha[pop[i]], 0, 0, Rmax[pop[i]], 0, 0,
                    S_sim[i], 0, 0, A[i], b) * exp(err_sim[i]);
  }
}
