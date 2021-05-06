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
  real pexp_lpdf(real y, real mu, real sigma_R, real shape) {
    return(-(fabs(y - mu)/sigma_R)^shape);
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
    if(cond[i]) which_cond = i;
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
  int<lower=0> N_X_alpha;              // number of intrinsic productivity covariates
  matrix[N,N_X_alpha] X_alpha;         // intrinsic productivity covariates
  int<lower=0> N_X_Rmax;               // number of maximum recruitment covariates
  matrix[N,N_X_Rmax] X_Rmax;           // maximum recruitment covariates
  int<lower=0> N_X_R;                  // number of recruitment covariates
  matrix[N,N_X_R] X_R;                 // brood-year productivity covariates
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
  int<lower=0> n_W_obs[N_H];           // count of wild spawners in samples
  int<lower=0> n_H_obs[N_H];           // count of hatchery spawners in samples
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years, not including fwd simulations
  int<lower=1,upper=N> N_year_all;         // total number of years, including forward simulations
  int<lower=1> pop_year_indx[N];           // index of years within each pop, starting at 1
  int<lower=2> ages[N_age];                // adult ages
  int<lower=1> min_age;                    // minimum adult age
  vector[N_age] ones_N_age = rep_vector(1,N_age); // for rowsums of p matrix 
  vector[N] ones_N = rep_vector(1,N);      // for elementwise inverse of rowsums 
  int<lower=0> n_HW_obs[N_H];              // total sample sizes for H/W frequencies
  int<lower=0,upper=N> fwd_init_indx[N_fwd,N_age]; // links "fitted" brood years to recruits in forward sims
  real mu_mu_Rmax = quantile(log(S_obs[which_S_obs]), 0.9); // prior mean of mu_Rmax
  real sigma_mu_Rmax = sd(log(S_obs[which_S_obs])); // prior SD of mu_Rmax
  vector[max_age*N_pop] mu_S_init;         // prior mean of total spawner abundance in years 1:max_age
  real sigma_S_init = 2*sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_age
  matrix[N_age,max_age*N_pop] mu_q_init;   // prior counts of wild spawner age distns in years 1:max_age
  
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
  
  for(i in 1:max_age)
  {
    int N_orphan_age = N_age - max(i - min_age, 0); // number of orphan age classes
    int N_amalg_age = N_age - N_orphan_age + 1;     // number of amalgamated age classes
    
    for(j in 1:N_pop)
    {
      int ii = (j - 1)*max_age + i; // index into S_init, q_init

      // S_init prior mean that scales observed log-mean by fraction of orphan age classes
      mu_S_init[ii] = mean(log(S_obs[which_S_obs])) + log(N_orphan_age) - log(N_age);
      
      // prior on q_init that implies q_orphan ~ Dir(1)
      mu_q_init[,ii] = append_row(rep_vector(1.0/N_amalg_age, N_amalg_age), 
                                  rep_vector(1, N_orphan_age - 1));
    }
  }
}

parameters {
  // recruitment
  real mu_alpha;                         // hyper-mean log intrinsic productivity
  vector[N_X_alpha] beta_alpha;          // regression coefs for log alpha
  real<lower=0> sigma_alpha;             // hyper-SD log intrinsic productivity
  vector[N_pop] zeta_alpha;              // log intrinsic prod (Z-scores)
  real mu_Rmax;                          // hyper-mean log maximum recruitment
  vector[N_X_Rmax] beta_Rmax;            // regression coefs for log Rmax
  real<lower=0> sigma_Rmax;              // hyper-SD log maximum recruitment
  vector[N_pop] zeta_Rmax;               // log maximum recruitment (Z-scores)
  real<lower=-1,upper=1> rho_alphaRmax;  // correlation between log(alpha) and log(Rmax)
  vector[N_X_R] beta_R;                  // regression coefs for log recruitment
  real<lower=-1,upper=1> rho_R;          // AR(1) coef for log productivity anomalies
  real<lower=0> sigma_year_R;            // hyper-SD of brood year log productivity anomalies
  vector[N_year_all] zeta_year_R;        // log brood year productivity anomalies (Z-scores)
  real<lower=0> sigma_R;                 // unique process error SD
  vector[N] zeta_R;                      // log true recruit abundance (not density) by brood year (z-scores)
  // spawner age structure 
  simplex[N_age] mu_p;                   // among-pop mean of age distributions
  vector<lower=0>[N_age-1] sigma_pop_p;  // among-pop SD of mean log-ratio age distributions
  cholesky_factor_corr[N_age-1] L_pop_p; // Cholesky factor of among-pop correlation matrix of mean log-ratio age distns
  matrix[N_pop,N_age-1] zeta_pop_p;      // population mean log-ratio age distributions (Z-scores)
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
  vector<lower=0>[N] alpha_Xbeta;        // intrinsic productivity including covariate effects
  vector<lower=0>[N_pop] Rmax;           // maximum recruitment 
  vector<lower=0>[N] Rmax_Xbeta;         // maximum recruitment including covariate effects
  vector[N_year_all] eta_year_R;         // log brood year productivity anomalies
  vector<lower=0>[N] R_hat;              // expected recruit abundance (not density) by brood year
  // vector<lower=0>[N] eta_Xbeta_epsilon_R; // recruitment process error including covariate effects
  vector<lower=0>[N] R;                  // true recruit abundance (not density) by brood year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                   // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;                // true total wild spawner abundance
  vector[N] S_H;                         // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                  // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  row_vector[N_age-1] mu_alr_p;          // mean of log-ratio cohort age distributions
  matrix[N_pop,N_age-1] mu_pop_alr_p;    // population mean log-ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p;    // cohort age distributions
  matrix<lower=0,upper=1>[N,N_age] q;    // true spawner age distributions
  
  // Multivariate Matt trick for [log(alpha), log(Rmax)]
  {
    matrix[2,2] L_alphaRmax;        // Cholesky factor of corr matrix of log(alpha), log(Rmax)
    matrix[N_pop,2] zeta_alphaRmax; // [log(alpha), log(Rmax)] random effects (z-scored)
    matrix[N_pop,2] eta_alphaRmax;  // [log(alpha), log(Rmax)] random effects
    vector[2] sigma_alphaRmax;      // SD vector of [log(alpha), log(Rmax)]
    
    L_alphaRmax[1,1] = 1;
    L_alphaRmax[2,1] = rho_alphaRmax;
    L_alphaRmax[1,2] = 0;
    L_alphaRmax[2,2] = sqrt(1 - rho_alphaRmax^2);
    sigma_alphaRmax[1] = sigma_alpha;
    sigma_alphaRmax[2] = sigma_Rmax;
    zeta_alphaRmax = append_col(zeta_alpha, zeta_Rmax);
    eta_alphaRmax = diag_pre_multiply(sigma_alphaRmax, L_alphaRmax * zeta_alphaRmax')';
    alpha = exp(mu_alpha + eta_alphaRmax[,1]);
    alpha_Xbeta = alpha[pop] .* exp(mat_lmult(X_alpha, beta_alpha));
    Rmax = exp(mu_Rmax + eta_alphaRmax[,2]);
    Rmax_Xbeta = Rmax[pop] .* exp(mat_lmult(X_Rmax, beta_Rmax));
  }
  
  // AR(1) model for eta_year_R
  eta_year_R[1] = zeta_year_R[1]*sigma_year_R/sqrt(1 - rho_R^2); // initial anomaly
  for(i in 2:N_year_all)
    eta_year_R[i] = rho_R*eta_year_R[i-1] + zeta_year_R[i]*sigma_year_R;
  // constrain "fitted" log anomalies to sum to 0
  eta_year_R = eta_year_R - mean(head(eta_year_R, N_year));

  // Total recruitment process error
  // eta_Xbeta_epsilon_R = exp(eta_year_R[year] + mat_lmult(X_R, beta_R) + sigma_R*zeta_R);
  
  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;
  
  // Multivariate Matt trick for age vectors (pop-specific mean and within-pop, time-varying)
  mu_alr_p = to_row_vector(log(head(mu_p, N_age-1)) - log(mu_p[N_age]));
  mu_pop_alr_p = rep_matrix(mu_alr_p,N_pop) + diag_pre_multiply(sigma_pop_p, L_pop_p * zeta_pop_p')';
  // Inverse log-ratio (softmax) transform of cohort age distn
  {
    matrix[N,N_age-1] alr_p = mu_pop_alr_p[pop,] + diag_pre_multiply(sigma_p, L_p * zeta_p')';
    matrix[N,N_age] exp_alr_p = append_col(exp(alr_p), ones_N);
    p = diag_pre_multiply(ones_N ./ (exp_alr_p * ones_N_age), exp_alr_p);
  }
  
  // Calculate true total wild and hatchery spawners and spawner age distribution
  // and predict recruitment from brood year i
  for(i in 1:N)
  {
    row_vector[N_age] S_W_a; // true wild spawners by age
    int ii; // index into S_init and q_init
    // number of orphan age classes <lower=0,upper=N_age>
    int N_orphan_age = max(N_age - max(pop_year_indx[i] - min_age, 0), N_age); 
    vector[N_orphan_age] q_orphan; // orphan age distribution (amalgamated simplex)

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
    R_hat[i] = SR(SR_fun, alpha_Xbeta[i], Rmax_Xbeta[i], S[i], A[i]);
    R[i] = R_hat[i] * exp(eta_year_R[year[i]] + dot_product(X_R[i], beta_R) + sigma_R*zeta_R[i]);
    // R[i] = R_hat[i] * eta_Xbeta_epsilon_R[i];
  }
}

model {
  vector[N_B] log_B_take; // log of true broodstock take when B_take_obs > 0
  
  // Priors
  
  // recruitment
  mu_alpha ~ normal(2,5);
  beta_alpha ~ normal(0,5);
  sigma_alpha ~ normal(0,3);
  mu_Rmax ~ normal(mu_mu_Rmax, sigma_mu_Rmax);
  beta_Rmax ~ normal(0,5);
  sigma_Rmax ~ normal(0,3);
  zeta_alpha ~ std_normal();  // [log(alpha), log(Rmax)] ~ MVN([mu_alpha, mu_Rmax], D*R_aRmax*D),
  zeta_Rmax ~ std_normal();   // where D = diag_matrix(sigma_alpha, sigma_Rmax)
  beta_R ~ normal(0,5);
  rho_R ~ pexp(0,0.85,50);    // mildly regularize to ensure stationarity
  zeta_year_R ~ std_normal(); // eta_year_R[i] ~ N(rho_R*eta_year_R[i-1], sigma_year_R)
  sigma_year_R ~ normal(0,3);
  sigma_R ~ normal(0,3);
  zeta_R ~ std_normal();      // total recruits: R ~ lognormal(log(R_hat), sigma_R)

  // spawner age structure
  for(i in 1:(N_age-1))
  {
    sigma_pop_p[i] ~ normal(0,2);
    sigma_p[i] ~ normal(0,2); 
  }
  L_pop_p ~ lkj_corr_cholesky(1);
  L_p ~ lkj_corr_cholesky(1);
  // pop mean age probs logistic MVN: 
  // mu_pop_alr_p[i,] ~ MVN(mu_alr_p,D*R_pop_p*D), where D = diag_matrix(sigma_pop_p)
  to_vector(zeta_pop_p) ~ std_normal();
  // age probs logistic MVN: 
  // alr_p[i,] ~ MVN(mu_pop_alr_p[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ std_normal();

  // removals
  log_B_take = log(S_W[which_B]) + log1m(q[which_B,1]) + logit(B_rate); // B_take = S_W*(1 - q[,1])*B_rate/(1 - B_rate)
  B_take_obs ~ lognormal(log_B_take, 0.05); // penalty to force pred and obs broodstock take to match 

  // initial spawners and wild spawner age distribution
  // (accounting for amalgamation of q_init to q_orphan)
  S_init ~ lognormal(mu_S_init, sigma_S_init);
  {
    matrix[N_age,max_age*N_pop] q_init_mat;
    
    for(j in 1:size(q_init)) q_init_mat[,j] = q_init[j];
    target += sum((mu_q_init - 1) .* log(q_init_mat)); // q_init[i] ~ Dir(mu_q_init[,i])
  }

  // spawner observation error
  tau ~ normal(0,1);
  
  // Observation model
  
  // total spawners (of observed ages)
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau); 
  n_H_obs ~ binomial(n_HW_obs, p_HOS); // counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));  // obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  // NOTE: DOES NOT INCLUDE COVARIATE EFFECTS //
  corr_matrix[N_age-1] R_pop_p;     // among-pop correlation matrix of mean log-ratio age distns 
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
  
  R_pop_p = multiply_lower_tri_self_transpose(L_pop_p);
  R_p = multiply_lower_tri_self_transpose(L_p);
  
  // Calculate true total wild and hatchery spawners and spawner age distribution
  // and simulate recruitment from brood year i
  // (Note that if N_fwd == 0, this block will not execute)
  for(i in 1:N_fwd)
  {
    vector[N_age-1] alr_p_fwd;     // alr(p_fwd[i,])'
    row_vector[N_age] S_W_a_fwd;   // true wild spawners by age
    
    // Inverse log-ratio transform of cohort age distn
    alr_p_fwd = multi_normal_cholesky_rng(to_vector(mu_pop_alr_p[pop_fwd[i],]), L_p);
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
    R_hat_fwd[i] = SR(SR_fun, alpha[pop_fwd[i]], Rmax[pop_fwd[i]], S_fwd[i], A_fwd[i]);
    R_fwd[i] = lognormal_rng(log(R_hat_fwd[i]) + eta_year_R[year_fwd[i]], sigma_R);
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
