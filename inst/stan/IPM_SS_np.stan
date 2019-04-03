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
}

data {
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  int<lower=1> N;                      // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         // population identifier
  int<lower=1,upper=N> year[N];        // brood year identifier
  int<lower=0> N_X;                    // number of productivity covariates
  matrix[max(year),N_X] X;             // brood-year productivity covariates (if none, use vector of zeros)
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs 
  int<lower=1,upper=N> which_S_obs[N_S_obs]; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed annual total spawner abundance (not density)
  int<lower=2> N_age;                  // number of adult age classes
  int<lower=2> max_age;                // maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  // observed wild spawner age frequencies (all zero row = NA)  
  int<lower=0,upper=N> N_H;            // number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[N_H];   // years with p_HOS > 0
  int<lower=0> n_W_obs[N_H];           // count of wild spawners in samples (assumes no NAs)
  int<lower=0> n_H_obs[N_H];           // count of hatchery spawners in samples (assumes no NAs)
  vector[N] A;                         // habitat area associated with each spawner abundance obs
  vector[N] F_rate;                    // fishing mortality rate of wild adults (no fishing on jacks)
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  int<lower=1,upper=N> which_B[N_B];   // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
}

transformed data {
  int<lower=1,upper=N> N_pop;     // number of populations
  int<lower=1,upper=N> N_year;    // number of years
  int<lower=2> ages[N_age];       // adult ages
  int<lower=1> pop_year_indx[N];  // index of years within each pop, starting at 1
  int<lower=0> n_HW_obs[N_H];     // total sample sizes for H/W frequencies
  
  N_pop = max(pop);
  N_year = max(year);
  for(a in 1:N_age)
    ages[a] = max_age - N_age + a;
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
}

parameters {
  vector<lower=0>[N_pop] alpha;           // intrinsic productivity
  vector<lower=0>[N_pop] Rmax;            // asymptotic recruitment
  matrix[N_pop,N_X] beta;                 // regression coefs for log productivity anomalies
  vector<lower=-1,upper=1>[N_pop] rho;    // AR(1) coefs for log productivity anomalies
  vector<lower=0>[N_pop] sigma;           // process error SDs
  simplex[N_age] mu_p[N_pop];             // population mean age distributions
  matrix<lower=0>[N_pop,N_age-1] sigma_p; // log-ratio cohort age distribution SDs
  cholesky_factor_corr[N_age-1] L_p[N_pop]; // Cholesky factors of correlation matrices of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;               // log-ratio cohort age distribution errors (Z-scores)
  vector<lower=0>[max_age*N_pop] S_init;  // true total spawner abundance in years 1:max_age
  simplex[N_age] q_init[max_age*N_pop];   // true wild spawner age distributions in years 1:max_age
  vector<lower=0,upper=1>[N_H] p_HOS;     // true p_HOS in years which_H
  vector[N] zeta_R;                       // recruitment process errors (z-scored)
  vector<lower=0,upper=1>[N_B] B_rate;    // true broodstock take rate when B_take > 0
  vector<lower=0>[N_pop] tau;             // observation error SDs of total spawners
}

transformed parameters {
  vector<lower=0>[N] S_W;             // true total wild spawner abundance
  vector[N] S_H;                      // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;               // true total spawner abundance
  matrix[N_pop,N_age-1] gamma;        // population mean log ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p; // cohort age distributions
  matrix<lower=0,upper=1>[N,N_age] q; // true spawner age distributions
  vector[N] p_HOS_all;                // true p_HOS in all years (can == 0)
  vector<lower=0>[N] R_hat;           // expected recruit abundance (not density) by brood year
  vector[N] epsilon_R;                // process error in recruit abundance by brood year 
  vector<lower=0>[N] R;               // true recruit abundance (not density) by brood year
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  
  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;
  
  // Log-ratio transform of pop-specific mean cohort age distributions
  for(j in 1:N_pop)
    gamma[j,] = to_row_vector(log(mu_p[j,1:(N_age-1)]) - log(mu_p[j,N_age]));
  
  // Calculate true total wild and hatchery spawners and spawner age distribution
  // and predict recruitment from brood year i
  for(i in 1:N)
  {
    row_vector[N_age] alr_p; // temp variable: alr(p[i,])
    row_vector[N_age] S_W_a; // temp variable: true wild spawners by age
    int ii;                  // temp variable: index into S_init and q_init
    
    // Multivariate Matt trick for within-pop, time-varying age vectors
    alr_p = rep_row_vector(0,N_age);
    alr_p[1:(N_age-1)] = gamma[pop[i],] + sigma_p[pop[i],] .* (L_p[pop[i]] * zeta_p[i,]')';
    alr_p = exp(alr_p);
    p[i,] = alr_p/sum(alr_p);
    
    if(pop_year_indx[i] <= max_age)
    {
      // use initial values
      ii = (pop[i] - 1)*max_age + pop_year_indx[i];
      S_W[i] = S_init[ii]*(1 - p_HOS_all[i]);        
      S_H[i] = S_init[ii]*p_HOS_all[i];
      q[i,] = to_row_vector(q_init[ii,]);
      S_W_a = S_W[i]*q[i,];
    }
    else
    {
      for(a in 1:N_age)
        S_W_a[a] = R[i-ages[a]]*p[i-ages[a],a];
      // catch and broodstock removal (assumes no take of age 1)
      S_W_a[2:N_age] = S_W_a[2:N_age]*(1 - F_rate[i])*(1 - B_rate_all[i]);
      S_W[i] = sum(S_W_a);
      S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
      q[i,] = S_W_a/S_W[i];
    }
    
    S[i] = S_W[i] + S_H[i];
    R_hat[i] = A[i] * SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S[i], A[i]);
    if(pop_year_indx[i] == 1) // initial process error
      epsilon_R[i] = zeta_R[i]*sigma[pop[i]]/sqrt(1 - rho[pop[i]]^2);
    else
      epsilon_R[i] = rho[pop[i]]*epsilon_R[i-1] + zeta_R[i]*sigma[pop[i]];
    R[i] = R_hat[i]*exp(dot_product(X[year[i],], beta[pop[i],]) + epsilon_R[i]);
  }
}

model {
  vector[N_B] B_take; // true broodstock take when B_take_obs > 0
  
  // Priors
  alpha ~ lognormal(2,2);
  Rmax ~ lognormal(2,3);
  to_vector(beta) ~ normal(0,5);
  for(j in 1:N_pop)
  {
    rho[j] ~ pexp(0,0.85,20);  // mildly regularize rho to ensure stationarity
    sigma[j] ~ pexp(0,1,10);
    tau[j] ~ pexp(1,0.85,30);  // rule out tau < 0.1 to avoid divergences 
    L_p[j] ~ lkj_corr_cholesky(3);
  }
  to_vector(sigma_p) ~ normal(0,5);
  S_init ~ lognormal(0,5);
  B_take = B_rate .* S_W[which_B] .* (1 - q[which_B,1]) ./ (1 - B_rate);
  B_take_obs ~ lognormal(log(B_take), 0.1); // penalty to force pred and obs broodstock take to match 

  // Hierarchical priors
  // age probs logistic MVN: alr_p[i,] ~ MVN(gamma[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ normal(0,1);
  
  // Process model
  zeta_R ~ normal(0,1); // total recruits: R ~ lognormal(log(R_hat), sigma)

  // Observation model
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau);  // observed total spawners
  n_H_obs ~ binomial(n_HW_obs, p_HOS); // obs counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));  // obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  corr_matrix[N_age-1] R_p[N_pop]; // correlation matrices of within-pop cohort log-ratio age distns
  vector[N] LL_S_obs;              // pointwise log-likelihood of total spawners
  vector[N_H] LL_n_H_obs;          // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_age_obs;          // pointwise log-likelihood of wild age frequencies
  vector[N] LL;                    // total pointwise log-likelihood                              
  
  for(j in 1:N_pop)
    R_p[j] = multiply_lower_tri_self_transpose(L_p[j]);
    
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
