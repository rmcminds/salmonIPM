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
  real pexp_lpdf(vector y, real mu, real sigma_R, real shape) {
    vector[num_elements(y)] LL;

    for(i in 1:num_elements(LL))
      LL[i] = -pow(fabs(y[i] - mu)/sigma_R, shape);

    return(sum(LL));
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
  vector[N] A;                         // habitat area associated with each spawner abundance obs
  // recruitment
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  int<lower=0> N_X;                    // number of productivity covariates
  matrix[max(year),N_X] X;             // brood-year productivity covariates
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
  // fishery and hatchery removals
  vector[N] F_rate;                    // fishing mortality rate of wild adults (no fishing on jacks)
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  int<lower=1,upper=N> which_B[N_B];   // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  int<lower=2> ages[N_age];                // adult ages
  int<lower=1> min_age;                    // minimum adult age
  int<lower=1> pop_year_indx[N];           // index of years within each pop, starting at 1
  int<lower=0> n_HW_obs[N_H];              // total sample sizes for H/W frequencies
  real mu_Rmax = quantile(log(S_obs[which_S_obs]), 0.9); // prior log-mean of Rmax
  real sigma_Rmax = sd(log(S_obs[which_S_obs])); // prior log-SD of Rmax
  vector[max_age*N_pop] mu_S_init;         // prior mean of total spawner abundance in years 1:max_age
  real sigma_S_init = 2*sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_age
  matrix[N_age,max_age*N_pop] mu_q_init;   // prior counts of wild spawner age distns in years 1:max_age

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
  vector<lower=0>[N_pop] alpha;           // intrinsic productivity
  vector<lower=0>[N_pop] Rmax;            // asymptotic recruitment
  matrix[N_pop,N_X] beta_R;               // regression coefs for log productivity anomalies
  vector<lower=-1,upper=1>[N_pop] rho_R;  // AR(1) coefs for log productivity anomalies
  vector<lower=0>[N_pop] sigma_R;         // process error SDs
  vector[N] zeta_R;                       // recruitment process errors (z-scored)
  // spawner age structure
  simplex[N_age] mu_p[N_pop];             // population mean age distributions
  matrix<lower=0>[N_pop,N_age-1] sigma_p; // log-ratio cohort age distribution SDs
  cholesky_factor_corr[N_age-1] L_p[N_pop]; // Cholesky factors of correlation matrices of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;               // log-ratio cohort age distribution errors (Z-scores)
  // H/W composition, removals
  vector<lower=0,upper=1>[N_H] p_HOS;     // true p_HOS in years which_H
  vector<lower=0,upper=1>[N_B] B_rate;    // true broodstock take rate when B_take > 0
  // initial spawners, observation error
  vector<lower=0>[max_age*N_pop] S_init;  // true total spawner abundance in years 1:max_age
  simplex[N_age] q_init[max_age*N_pop];   // true wild spawner age distributions in years 1:max_age
  vector<lower=0>[N_pop] tau;             // observation error SDs of total spawners
}

transformed parameters {
  // recruitment
  vector<lower=0>[N] R_hat;            // expected recruit abundance (not density) by brood year
  vector[N] epsilon_R;                 // process error in recruit abundance by brood year
  vector<lower=0>[N] R;                // true recruit abundance (not density) by brood year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                 // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;              // true total wild spawner abundance
  vector[N] S_H;                       // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  vector[N_age-1] mu_alr_p[N_pop];     // population mean log ratio age distributions
  simplex[N_age] p[N];                 // cohort age distributions
  matrix<lower=0,upper=1>[N,N_age] q;  // true spawner age distributions

  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;

  // Log-ratio transform of pop-specific mean cohort age distributions
  for(j in 1:N_pop)
    mu_alr_p[j] = log(head(mu_p[j], N_age-1)) - log(tail(mu_p[j], 1));

  // Calculate true total wild and hatchery spawners and spawner age distribution
  // and predict recruitment from brood year i
  for(i in 1:N)
  {
    vector[N_age] alr_p;     // alr(p[i,])
    row_vector[N_age] S_W_a; // true wild spawners by age
    int ii;                  // index into S_init and q_init
    // number of orphan age classes <lower=0,upper=N_age>
    int N_orphan_age = max(N_age - max(pop_year_indx[i] - min_age, 0), N_age); 
    vector[N_orphan_age] q_orphan; // orphan age distribution (amalgamated simplex)

    // Multivariate Matt trick for within-pop, time-varying age vectors
    alr_p = rep_vector(0,N_age);
    alr_p[1:(N_age-1)] = mu_alr_p[pop[i]] + sigma_p[pop[i],]' .* (L_p[pop[i]] * zeta_p[i,]');
    alr_p = exp(alr_p);
    p[i] = alr_p / sum(alr_p);
    
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
    R_hat[i] = SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S[i], A[i]);
    if(pop_year_indx[i] == 1) // initial process error
      epsilon_R[i] = zeta_R[i]*sigma_R[pop[i]]/sqrt(1 - rho_R[pop[i]]^2);
    else
      epsilon_R[i] = rho_R[pop[i]]*epsilon_R[i-1] + zeta_R[i]*sigma_R[pop[i]];
    R[i] = R_hat[i]*exp(dot_product(X[year[i],], beta_R[pop[i],]) + epsilon_R[i]);
  }
}

model {
  vector[N_B] log_B_take; // log of true broodstock take when B_take_obs > 0
  
  // Priors
  
  // recruitment
  alpha ~ lognormal(2.0,2.0);
  Rmax ~ lognormal(mu_Rmax, sigma_Rmax);
  to_vector(beta_R) ~ normal(0,5);
  rho_R ~ pexp(0,0.85,20); // mildly regularize rho_R to ensure stationarity
  sigma_R ~ normal(0,5);
  zeta_R ~ std_normal();   // total recruits: R ~ lognormal(log(R_hat), sigma_R)

  // spawner age structure
  to_vector(sigma_p) ~ normal(0,5);
  for(j in 1:N_pop)
    L_p[j] ~ lkj_corr_cholesky(3);
  // age probs logistic MVN: alr_p[i,] ~ MVN(mu_alr_p[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
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
  tau ~ pexp(1,0.85,30);  // rule out tau < 0.1 to avoid divergences

  // Observation model
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau[pop[which_S_obs]]);  // observed total spawners
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
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau[pop[which_S_obs]]);
  LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
  LL_n_H_obs = rep_vector(0,N_H);
  for(i in 1:N_H)
    LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  LL = LL_S_obs + LL_n_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
