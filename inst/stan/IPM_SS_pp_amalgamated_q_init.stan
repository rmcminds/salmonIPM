functions {
  // spawner-recruit functions
  real SR(int SR_fun, real alpha, real Rmax, real S, real A) {
    real R;
    
    if(SR_fun == 1)      // discrete exponential
    R = alpha*S/A;
    else if(SR_fun == 2) // Beverton-Holt
    R = alpha*S/(A + alpha*S/Rmax);
    else if(SR_fun == 3) // Ricker
    R = alpha*(S/A)*exp(-alpha*S/(A*e()*Rmax));
    
    return(R);
  }
  
  // Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(real y, real mu, real sigma, real shape) {
    return(-(fabs(y - mu)/sigma)^shape);
  }
  
  // convert matrix to array of column vectors
  vector[] matrix_to_array(matrix m) {
    vector[2] arr[cols(m)];
    
    for(i in 1:cols(m))
      arr[i] = col(m,i);
    return(arr);
  }
  
  // Vectorized logical equality
  int[] veq(int[] x, int y) {
    int xeqy[size(x)];
    for(i in 1:size(x))
      xeqy[i] = x[i] == y;
    return(xeqy);
  }
  
  // Vectorized logical &&
    int[] vand(int[] cond1, int[] cond2) {
      int cond1_and_cond2[size(cond1)];
      for(i in 1:size(cond1))
        cond1_and_cond2[i] = cond1[i] && cond2[i];
      return(cond1_and_cond2);
    }
  
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
  
  // Equivalent of R: which(cond), where sum(cond) == 1
  int which(int[] cond) {
    int which_cond;
    for(i in 1:size(cond))
      if(cond[i])
        which_cond = i;
      return(which_cond);
  }
  
  // Left multiply vector by matrix
  // works even if size is zero
  vector mat_lmult(matrix X, vector v)
  {
    vector[rows(X)] Xv;
    Xv = rows_dot_product(X, rep_matrix(to_row_vector(v), rows(X)));
    return(Xv); 
  }
}

data {
  // info for observed data
  int<lower=1> N;                      // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         // population identifier
  int<lower=1,upper=N> year[N];        // brood year identifier
  // info for forward simulations
  int<lower=0> N_fwd;                  // total number of cases in forward simulations
  int<lower=1,upper=N> pop_fwd[N_fwd]; // population identifier for forward simulations
  int<lower=1,upper=N+N_fwd> year_fwd[N_fwd]; // brood year identifier for forward simulations
  vector<lower=0>[N_fwd] A_fwd;        // habitat area for each forward simulation
  vector<lower=0,upper=1>[N_fwd] F_rate_fwd; // fishing mortality for forward simulations
  vector<lower=0,upper=1>[N_fwd] B_rate_fwd; // broodstock take rate for forward simulations
  vector<lower=0,upper=1>[N_fwd] p_HOS_fwd; // p_HOS for forward simulations
  // recruitment
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  vector<lower=0>[N] A;                // habitat area associated with each spawner abundance obs
  int<lower=0> N_X;                    // number of productivity covariates
  matrix[max(append_array(year,year_fwd)),N_X] X; // brood-year productivity covariates
  // fishery and hatchery removals
  vector<lower=0,upper=1>[N] F_rate;   // fishing mortality of wild adults
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  int<lower=1,upper=N> which_B[N_B];   // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
  // spawner abundance
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs 
  int<lower=1,upper=N> which_S_obs[N_S_obs]; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed annual total spawner abundance (not density)
  // spawner age structure
  int<lower=2> N_age;                  // number of adult age classes
  int<lower=2> max_age;                // maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  // observed wild spawner age frequencies (all zero row = NA)  
  // H/W composition
  int<lower=0,upper=N> N_H;            // number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[N_H];   // years with p_HOS > 0
  int<lower=0> n_W_obs[N_H];           // count of wild spawners in samples (assumes no NAs)
  int<lower=0> n_H_obs[N_H];           // count of hatchery spawners in samples (assumes no NAs)
}

transformed data {
  int<lower=1,upper=N> N_pop;        // number of populations
  int<lower=1,upper=N> N_year;       // number of years, not including forward simulations
  int<lower=1,upper=N> N_year_all;   // total number of years, including forward simulations
  int<lower=2> ages[N_age];          // adult ages
  int<lower=1> min_age;              // minimum adult age
  int<lower=0> n_HW_obs[N_H];        // total sample sizes for H/W frequencies
  int<lower=1> pop_year_indx[N];     // index of years within each pop, starting at 1
  int<lower=0,upper=N> fwd_init_indx[N_fwd,N_age]; // links "fitted" brood years to recruits in forward sims
  
  N_pop = max(pop);
  N_year = max(year);
  N_year_all = max(append_array(year, year_fwd));
  for(a in 1:N_age)
    ages[a] = max_age - N_age + a;
  min_age = min(ages);  
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
  
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  
  fwd_init_indx = rep_array(0, N_fwd, N_age);
  for(i in 1:N_fwd)
  {
    for(a in 1:N_age)
    {
      if(year_fwd[i] - ages[a] < min(rsub(year_fwd, veq(pop_fwd, pop_fwd[i]))))
        fwd_init_indx[i,a] = which(vand(veq(pop, pop_fwd[i]), veq(year, year_fwd[i] - ages[a])));
    }
  }
}

parameters {
  // recruitment
  real mu_alpha;                         // hyper-mean log intrinsic productivity
  real<lower=0> sigma_alpha;             // hyper-SD log intrinsic productivity
  vector[N_pop] zeta_alpha;              // log intrinsic prod (Z-scores)
  real mu_Rmax;                          // hyper-mean log asymptotic recruitment
  real<lower=0> sigma_Rmax;              // hyper-SD log asymptotic recruitment
  vector[N_pop] zeta_Rmax;               // log asymptotic recruitment (Z-scores)
  real<lower=-1,upper=1> rho_alphaRmax;  // correlation between log(alpha) and log(Rmax)
  vector[N_X] beta_phi;                  // regression coefs for log productivity anomalies
  real<lower=-1,upper=1> rho_phi;        // AR(1) coef for log productivity anomalies
  real<lower=0> sigma_phi;               // hyper-SD of brood year log productivity anomalies
  vector[N_year_all] zeta_phi;           // log brood year productivity anomalies (Z-scores)
  real<lower=0> sigma;                   // unique process error SD
  vector[N] zeta_R;                      // log true recruit abundance (not density) by brood year (z-scores)
  // spawner age structure 
  simplex[N_age] mu_p;                   // among-pop mean of age distributions
  vector<lower=0>[N_age-1] sigma_gamma;  // among-pop SD of mean log-ratio age distributions
  cholesky_factor_corr[N_age-1] L_gamma; // Cholesky factor of among-pop correlation matrix of mean log-ratio age distns
  matrix[N_pop,N_age-1] zeta_gamma;      // population mean log-ratio age distributions (Z-scores)
  vector<lower=0>[N_age-1] sigma_p;      // SD of log-ratio cohort age distributions
  cholesky_factor_corr[N_age-1] L_p;     // Cholesky factor of correlation matrix of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;              // log-ratio cohort age distributions (Z-scores)
  // H/W composition, removals
  vector<lower=0,upper=1>[N_H] p_HOS;    // true p_HOS in years which_H
  vector<lower=0,upper=1>[N_B] B_rate;   // true broodstock take rate when B_take > 0
  // initial spawners, observation error
  vector<lower=0>[max_age*N_pop] S_init; // true total spawner abundance in years 1-max_age
  simplex[N_age] q_init[max_age*N_pop];  // true wild spawner age distributions in years 1-max_age
  real<lower=0> tau;                     // observation error SD of total spawners
}

transformed parameters {
  // recruitment
  vector<lower=0>[N_pop] alpha;          // intrinsic productivity 
  vector<lower=0>[N_pop] Rmax;           // asymptotic recruitment 
  vector[N_year_all] phi;                // log brood year productivity anomalies
  vector<lower=0>[N] R_hat;              // expected recruit abundance (not density) by brood year
  vector<lower=0>[N] R;                  // true recruit abundance (not density) by brood year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                   // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;                // true total wild spawner abundance
  vector[N] S_H;                         // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                  // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  row_vector[N_age-1] mu_gamma;          // mean of log-ratio cohort age distributions
  matrix[N_pop,N_age-1] gamma;           // population mean log-ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p;    // cohort age distributions
  matrix<lower=0,upper=1>[N,N_age] q;    // true spawner age distributions
  
  // Multivariate Matt trick for [log(alpha), log(Rmax)]
  {
    matrix[2,2] L_alphaRmax;           // Cholesky factor of corr matrix of log(alpha), log(Rmax)
    matrix[N_pop,2] zeta_alphaRmax;    // [log(alpha), log(Rmax)] random effects (z-scored)
    matrix[N_pop,2] epsilon_alphaRmax; // [log(alpha), log(Rmax)] random effects
    vector[2] sigma_alphaRmax;         // SD vector of [log(alpha), log(Rmax)]
    
    L_alphaRmax[1,1] = 1;
    L_alphaRmax[2,1] = rho_alphaRmax;
    L_alphaRmax[1,2] = 0;
    L_alphaRmax[2,2] = sqrt(1 - rho_alphaRmax^2);
    sigma_alphaRmax[1] = sigma_alpha;
    sigma_alphaRmax[2] = sigma_Rmax;
    zeta_alphaRmax = append_col(zeta_alpha, zeta_Rmax);
    epsilon_alphaRmax = (diag_matrix(sigma_alphaRmax) * L_alphaRmax * zeta_alphaRmax')';
                         alpha = exp(mu_alpha + col(epsilon_alphaRmax,1));
                         Rmax = exp(mu_Rmax + col(epsilon_alphaRmax,2));
  }
  
  // AR(1) model for phi
  phi[1] = zeta_phi[1]*sigma_phi/sqrt(1 - rho_phi^2); // initial anomaly
  for(i in 2:N_year_all)
    phi[i] = rho_phi*phi[i-1] + zeta_phi[i]*sigma_phi;
  // constrain "fitted" log anomalies to sum to 0 (X should be centered)
  phi = phi - mean(phi[1:N_year]) + mat_lmult(X,beta_phi);
  
  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;
  
  // Multivariate Matt trick for age vectors (pop-specific mean and within-pop, time-varying)
  mu_gamma = to_row_vector(log(mu_p[1:(N_age-1)]) - log(mu_p[N_age]));
  gamma = rep_matrix(mu_gamma,N_pop) + (diag_matrix(sigma_gamma) * L_gamma * zeta_gamma')';
  p = append_col(gamma[pop,] + (diag_matrix(sigma_p) * L_p * zeta_p')', rep_vector(0,N));
  
  // Calculate true total wild and hatchery spawners and spawner age distribution
  // and predict recruitment from brood year i
  for(i in 1:N)
  {
    row_vector[N_age] exp_p; // exp(p[i,])
    row_vector[N_age] S_W_a; // true wild spawners by age
    int ii; // index into S_init and q_init
    // number of orphan age classes <lower=0,upper=N_age>
    int N_orphan_age = max(N_age - max(pop_year_indx[i] - min_age, 0), N_age); 
    vector[N_orphan_age] q_orphan; // orphan age distribution (amalgamated simplex)
    
    // Inverse log-ratio transform of cohort age distn
    // (built-in softmax function doesn't accept row vectors)
    exp_p = exp(p[i,]);
    p[i,] = exp_p/sum(exp_p);
    
    // Use initial values for orphan age classes, otherwise use process model
    if(pop_year_indx[i] <= max_age)
    {
      ii = (pop[i] - 1)*max_age + pop_year_indx[i];
      q_orphan = append_row(sum(head(q_init[ii], N_age - N_orphan_age + 1)), 
                            tail(q_init[ii], N_orphan_age - 1));
    }
    
    for(a in 1:N_age)
    {
      if(ages[a] < pop_year_indx[i])
        // Use recruitment process model
        S_W_a[a] = R[i-ages[a]]*p[i-ages[a],a];
      else
        // Use initial values
        S_W_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_orphan[a - (N_age - N_orphan_age)];
    }
    
    // catch and broodstock removal (assumes no take of age 1)
    S_W_a[2:N_age] = S_W_a[2:N_age]*(1 - F_rate[i])*(1 - B_rate_all[i]);
    S_W[i] = sum(S_W_a);
    S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
    S[i] = S_W[i] + S_H[i];
    q[i,] = S_W_a/S_W[i];
    
    // Recruitment
    R_hat[i] = A[i]*SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S[i], A[i]);
    R[i] = R_hat[i]*exp(phi[year[i]] + sigma*zeta_R[i]);
  }
}

model {
  vector[N_B] B_take; // true broodstock take when B_take_obs > 0
  
  // Priors
  
  // recruitment
  mu_alpha ~ normal(2,5);
  sigma_alpha ~ pexp(0,3,10);
  mu_Rmax ~ normal(0,10);
  sigma_Rmax ~ pexp(0,3,10);
  zeta_alpha ~ std_normal(); // [log(alpha), log(Rmax)] ~ MVN([mu_alpha, mu_Rmax], D*R_aRmax*D),
  zeta_Rmax ~ std_normal();  // where D = diag_matrix(sigma_alpha, sigma_Rmax)
  beta_phi ~ normal(0,5);
  rho_phi ~ pexp(0,0.85,50); // mildly regularize to ensure stationarity
  zeta_phi ~ std_normal();   // phi[i] ~ N(rho_phi*phi[i-1], sigma_phi)
  sigma_phi ~ pexp(0,2,10);
  sigma ~ pexp(0,1,10);
  zeta_R ~ std_normal();     // total recruits: R ~ lognormal(log(R_hat), sigma)

  // spawner age structure
  for(i in 1:(N_age-1))
  {
    sigma_gamma[i] ~ pexp(0,2,5);
    sigma_p[i] ~ pexp(0,2,5); 
  }
  L_gamma ~ lkj_corr_cholesky(1);
  L_p ~ lkj_corr_cholesky(1);
  // pop mean age probs logistic MVN: 
  // gamma[i,] ~ MVN(mu_gamma,D*R_gamma*D), where D = diag_matrix(sigma_gamma)
  to_vector(zeta_gamma) ~ std_normal();
  // age probs logistic MVN: 
  // alr_p[i,] ~ MVN(gamma[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ std_normal();

  // removals
  B_take = B_rate .* S_W[which_B] .* (1 - q[which_B,1]) ./ (1 - B_rate);
  B_take_obs ~ lognormal(log(B_take), 0.1); // penalty to force pred and obs broodstock take to match 

  // initial spawners and wild spawner age distribution
  // (accounting for amalgamation of q_init to q_orphan)
  for(i in 1:N)
  {
    if(pop_year_indx[i] <= max_age)
    {
      int N_orphan_age = N_age - max(pop_year_indx[i] - min_age, 0); // # orphan age classes
      int N_amalg_age = N_age - N_orphan_age + 1; // # amalgamated age classes
      int ii = (pop[i] - 1)*max_age + pop_year_indx[i]; // index into q_init
      
      S_init[ii] ~ lognormal(log(1.0*N_orphan_age/N_age), 10.0);
      
      // prior on q_init that implies q_orphan ~ Dir(1)
      q_init[ii] ~ dirichlet(append_row(rep_vector(1.0/N_amalg_age, N_amalg_age),
                                        rep_vector(1, N_orphan_age - 1)));
    }
  }

  // spawner observation error
  tau ~ pexp(0,1,10);
  
  // Observation model
  
  // total spawners (of observed ages)
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau); 
  n_H_obs ~ binomial(n_HW_obs, p_HOS); // counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));  // obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  corr_matrix[N_age-1] R_gamma;     // among-pop correlation matrix of mean log-ratio age distns 
  corr_matrix[N_age-1] R_p;         // correlation matrix of within-pop cohort log-ratio age distns 
  vector<lower=0>[N_fwd] S_W_fwd;   // true total wild spawner abundance in forward simulations
  vector[N_fwd] S_H_fwd;            // true total hatchery spawner abundance in forward simulations
  vector<lower=0>[N_fwd] S_fwd;     // true total spawner abundance in forward simulations
  matrix<lower=0,upper=1>[N_fwd,N_age] p_fwd; // cohort age distributions in forward simulations
  matrix<lower=0,upper=1>[N_fwd,N_age] q_fwd; // spawner age distributions in forward simulations
  vector<lower=0>[N_fwd] R_hat_fwd; // expected recruit abundance by brood year in forward simulations
  vector<lower=0>[N_fwd] R_fwd;     // true recruit abundance by brood year in forward simulations
  vector[N] LL_S_obs;               // pointwise log-likelihood of total spawners
  vector[N_H] LL_n_H_obs;           // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_age_obs;           // pointwise log-likelihood of wild age frequencies
  vector[N] LL;                     // total pointwise log-likelihood                              
  
  R_gamma = multiply_lower_tri_self_transpose(L_gamma);
  R_p = multiply_lower_tri_self_transpose(L_p);
  
  // Calculate true total wild and hatchery spawners and spawner age distribution
  // and simulate recruitment from brood year i
  // (Note that if N_fwd == 0, this block will not execute)
  for(i in 1:N_fwd)
  {
    vector[N_age-1] alr_p_fwd;     // alr(p_fwd[i,])'
    row_vector[N_age] S_W_a_fwd;   // true wild spawners by age

// Inverse log-ratio transform of cohort age distn
alr_p_fwd = multi_normal_cholesky_rng(to_vector(gamma[pop_fwd[i],]), L_p);
p_fwd[i,] = to_row_vector(softmax(append_row(alr_p_fwd,0)));

for(a in 1:N_age)
{
  if(fwd_init_indx[i,a] != 0)
  {
    // Use estimated values from previous cohorts
    S_W_a_fwd[a] = R[fwd_init_indx[i,a]]*p[fwd_init_indx[i,a],a];
  }
  else
  {
    S_W_a_fwd[a] = R_fwd[i-ages[a]]*p_fwd[i-ages[a],a];
  }
}

for(a in 2:N_age)  // catch and broodstock removal (assumes no take of age 1)
S_W_a_fwd[a] = S_W_a_fwd[a]*(1 - F_rate_fwd[i])*(1 - B_rate_fwd[i]);

S_W_fwd[i] = sum(S_W_a_fwd);
S_H_fwd[i] = S_W_fwd[i]*p_HOS_fwd[i]/(1 - p_HOS_fwd[i]);
q_fwd[i,] = S_W_a_fwd/S_W_fwd[i];
S_fwd[i] = S_W_fwd[i] + S_H_fwd[i];
R_hat_fwd[i] = A_fwd[i] * SR(SR_fun, alpha[pop_fwd[i]], Rmax[pop_fwd[i]], S_fwd[i], A_fwd[i]);
R_fwd[i] = lognormal_rng(log(R_hat_fwd[i]) + phi[year_fwd[i]], sigma);
                                                       }

LL_S_obs = rep_vector(0,N);
for(i in 1:N_S_obs)
  LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau); 
LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
LL_n_H_obs = rep_vector(0,N_H);
for(i in 1:N_H)
  LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
LL = LL_S_obs + LL_n_age_obs;
LL[which_H] = LL[which_H] + LL_n_H_obs;
}
