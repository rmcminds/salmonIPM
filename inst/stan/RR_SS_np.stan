functions {
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
  
  // Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(real y, real mu, real sigma, real shape) {
    return(-(fabs(y - mu)/sigma)^shape);
  }
    
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
}

data {
  int<lower=1> SR_fun;          // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  int<lower=1> N;               // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];  // population identifier
  int<lower=1,upper=N> year[N]; // brood year identifier
  int<lower=1,upper=N> N_fit;   // number of cases used in fitting (non-missing S and R)
  int<lower=1,upper=N> which_fit[N_fit]; // cases used in fitting
  vector<lower=0>[N] S;         // observed annual total spawner abundance (not density)
  vector<lower=0>[N] R;         // total natural recruit abundance (not density), including harvest and broodstock removals
  vector[N] A;                  // habitat area associated with each spawner abundance obs
  int<lower=0,upper=1> S_NA[N]; // logical indicating whether S is missing and should be simulated
  int<lower=0,upper=1> R_NA[N]; // logical indicating whether R is missing and should be simulated
  int<lower=2> N_age;           // number of adult age classes
  int<lower=2> max_age;         // maximum adult age
  matrix<lower=0,upper=1>[max(pop),N_age] p;  // average recruit age distributions for each pop 
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years, not incl fwd simulations
  int<lower=2> ages[N_age];                // adult ages
  real mu_mu_Rmax = quantile(log(R[which_fit]), 0.9);  // prior mean of mu_Rmax
  real sigma_mu_Rmax = 2*sd(log(R[which_fit]));  // prior SD of mu_Rmax
  
  for(a in 1:N_age)
    ages[a] = max_age - N_age + a;
}

parameters {
  vector<lower=0>[N_pop] alpha;         // intrinsic productivity
  vector<lower=0>[N_pop] Rmax;          // asymptotic recruitment
  vector<lower=-1,upper=1>[N_pop] rho;  // AR(1) coefs of residuals
  vector<lower=0>[N_pop] sigma;         // residual error SD
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
    R_hat[which_fit[i]] = SR(SR_fun, alpha[pop[which_fit[i]]], Rmax[pop[which_fit[i]]], 
                             S[which_fit[i]], A[which_fit[i]]);
    if(i==1 || pop[which_fit[i-1]] != pop[which_fit[i]])
    {
      R_ar1[which_fit[i]] = R_hat[which_fit[i]];
      sigma_ar1[which_fit[i]] = sigma[pop[which_fit[i]]]/sqrt(1 - rho[pop[which_fit[i]]]^2);
    }
    else
    {
      real err;    // temp variable: residual at the last non-missing observation
      int dt;      // temp variable: number of years since last non-missing observation
      real rho2j;  // temp variable: sum of powers of rho
      
      err = log(R[which_fit[i-1]]) - log(R_hat[which_fit[i-1]]);
      dt = which_fit[i] - which_fit[i-1];
      R_ar1[which_fit[i]] = R_hat[which_fit[i]]*exp(rho[pop[which_fit[i]]]^dt * err);
      
      rho2j = 0;
      for(j in 0:(dt-1))
        rho2j = rho2j + rho[pop[which_fit[i]]] ^ (2*j);
      sigma_ar1[which_fit[i]] = sigma[pop[which_fit[i]]] * sqrt(rho2j);
    }
  }
}

model {
  // Priors
  alpha ~ lognormal(2.0,2.0);
  Rmax ~ lognormal(mu_mu_Rmax, sigma_mu_Rmax);
  for(i in 1:N_pop)
  {
    rho[i] ~ pexp(0,0.85,20);   // mildly regularize rho to ensure stationarity
    sigma[i] ~ normal(0,3);
  }

  // Likelihood
  R[which_fit] ~ lognormal(log(R_ar1[which_fit]), sigma_ar1[which_fit]);
}

generated quantities {
  vector[N] S_sim;    // simulated spawners
  vector[N] R_sim;    // simulated recruits
  vector[N] err_sim;  // simulated AR(1) residual errors
  
  S_sim = S;
  R_sim = R;
  err_sim = rep_vector(0,N);
  
  for(i in 1:N)
  {
    if(S_NA[i] == 1)
    {
      if(i >= max_age && pop[i-max_age] == pop[i])
      {
        S_sim[i] = 0;
        for(a in 1:N_age)
          S_sim[i] = S_sim[i] + R_sim[i-ages[a]]*p[pop[i],a];
      }
    }
    
    if(i == 1 || pop[i-1] != pop[i])
      err_sim[i] = normal_rng(0, sigma[pop[i]]);
    else
      err_sim[i] = normal_rng(rho[pop[i]]*err_sim[i-1], sigma[pop[i]]);
    
    if(R_NA[i] == 1)
      R_sim[i] = SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S_sim[i], A[i])*exp(err_sim[i]);
  }
}
