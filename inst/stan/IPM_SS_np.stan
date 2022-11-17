functions {
  #include /include/SR.stan
  #include /include/gnormal_lpdf_vec.stan
  #include /include/quantile.stan
}

data {
  // info for observed data
  int<lower=1> N;                      // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         // population index
  int<lower=1,upper=N> year[N];        // brood year index
  // recruitment
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  vector<lower=0>[N] A;                // habitat area associated with each spawner abundance obs
  int<lower=0> K_alpha;                // number of intrinsic productivity covariates
  matrix[N,K_alpha] X_alpha;           // intrinsic productivity covariates
  real prior_alpha[2];                 // prior meanlog, sdlog for intrinsic productivity
  int<lower=0> K_Rmax;                 // number of maximum recruitment covariates
  matrix[N,K_Rmax] X_Rmax;             // maximum recruitment covariates
  real prior_Rmax[2];                  // prior meanlog, sdlog for maximum recruitment
  int<lower=0> K_R;                    // number of recruitment covariates
  row_vector[K_R] X_R[N];              // brood-year productivity covariates
  // kelt survival
  int<lower=0,upper=1> iter;           // is life cycle semelparous (0) or iteroparous (1)?
  int<lower=0> K_SS;                   // number of kelt survival covariates
  row_vector[K_SS] X_SS[N];            // kelt survival covariates
  real prior_mu_SS[2];                 // prior a, b for mean kelt survival
  // spawner abundance
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs
  int<lower=1,upper=N> which_S_obs[N_S_obs]; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed total spawner abundance (not density)
  real prior_tau[3];                   // prior mean, scale, shape for spawner observation error SD
  // spawner age structure
  int<lower=2> N_age;                  // number of (maiden) adult age classes
  int<lower=2> max_age;                // maximum (maiden) adult age
  matrix<lower=0>[N,iter ? N_age*2 : N_age] n_age_obs; // wild spawner [maiden | kelt] age frequencies
  vector<lower=0>[N_age] prior_mu_p;   // prior concentration for mean age distribution
  // H/W composition
  int<lower=0,upper=N> N_H;            // number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[N_H];   // years with p_HOS > 0
  int<lower=0> n_W_obs[N_H];           // count of wild spawners in samples
  int<lower=0> n_H_obs[N_H];           // count of hatchery spawners in samples
  // fishery and hatchery removals
  vector[N] F_rate;                    // fishing mortality rate of wild adults
  vector<lower=0,upper=1>[N_age + iter] age_F; // is age a (non)selected (0/1) by fishery?
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  int<lower=1,upper=N> which_B[N_B];   // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
  vector<lower=0,upper=1>[N_age + iter] age_B; // is age a (non)selected (0/1) in broodstock?
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  int<lower=2> ages[N_age];                // (maiden) adult ages
  int<lower=1> min_age;                    // minimum adult age
  int<lower=1> pop_year[N];                // index of years within each pop, starting at 1
  int<lower=0> n_HW_obs[N_H];              // total sample sizes for H/W frequencies
  vector[max_age*N_pop] mu_S_init;         // prior mean of total spawner abundance in years 1:max_age
  real sigma_S_init = sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_age
  matrix[N_age,max_age*N_pop] mu_q_init;   // prior counts of maiden age distributions in years 1:max_age
  
  for(a in 1:N_age) ages[a] = max_age - N_age + a;
  min_age = min(ages);  
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
  
  pop_year[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i]) pop_year[i] = 1;
    else pop_year[i] = pop_year[i-1] + 1;
  }
  
  for(i in 1:max_age)
  {
    int N_orphan_age = N_age - max(i - min_age, 0); // number of orphan (maiden) age classes
    int N_amalg_age = N_age - N_orphan_age + 1;     // number of amalgamated (maiden) age classes
    
    for(j in 1:N_pop)
    {
      int ii = (j - 1)*max_age + i; // index into S_init, q_init
      
      // S_init prior mean scales observed log-mean by fraction of orphan maiden age classes
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
  matrix[N_pop,K_alpha] beta_alpha;       // regression coefs for log alpha
  vector<lower=0>[N_pop] Rmax;            // maximum recruitment
  matrix[N_pop,K_Rmax] beta_Rmax;         // regression coefs for log Rmax
  matrix[N_pop,K_R] beta_R;               // regression coefs for log recruitment
  vector<lower=-1,upper=1>[N_pop] rho_R;  // AR(1) coefs for log productivity anomalies
  vector<lower=0>[N_pop] sigma_R;         // process error SDs
  vector[N] zeta_R;                       // recruitment process errors (Z-scored)
  // (maiden) spawner age structure
  simplex[N_age] mu_p[N_pop];             // population mean age distributions
  matrix<lower=0>[N_pop,N_age-1] sigma_p; // log-ratio cohort age distribution SDs
  cholesky_factor_corr[N_age-1] L_p[N_pop]; // Cholesky factors of correlation matrices of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;               // log-ratio cohort age distribution errors (Z-scores)
  // kelt survival
  vector<lower=0,upper=1>[iter*N_pop] mu_SS; // mean kelt (spawner-to-spawner) survival
  matrix[N_pop,K_SS] beta_SS;             // regression coefs for kelt survival
  vector<lower=-1,upper=1>[iter*N_pop] rho_SS; // AR(1) coefs for kelt survival
  vector<lower=0>[iter*N_pop] sigma_SS;   // kelt survival process error SDs
  vector[iter*N] zeta_SS;                 // kelt survival process errors (Z-scored)
  // H/W composition, removals
  vector<lower=0,upper=1>[N_H] p_HOS;     // true p_HOS in years which_H
  vector<lower=0,upper=1>[N_B] B_rate;    // true broodstock take rate when B_take > 0
  // initial spawners, observation error
  vector<lower=0>[max_age*N_pop] S_init;  // true total spawner abundance in years 1:max_age
  simplex[N_age] q_init[max_age*N_pop];   // true wild (maiden) age distribution in years 1:max_age
  simplex[N_age*2] q_iter_init[iter*N_pop]; // true wild [maiden | kelt] age distribution in year 1
  vector<lower=0>[N_pop] tau;             // observation error SDs of total spawners
}

transformed parameters {
  // recruitment
  vector<lower=0>[N] alpha_Xbeta;      // intrinsic productivity including covariate effects
  vector<lower=0>[N] Rmax_Xbeta;       // maximum recruitment including covariate effects
  vector<lower=0>[N] R_hat;            // expected recruit abundance (not density) by brood year
  vector[N] epsilon_R;                 // process error in recruit abundance by brood year
  vector<lower=0>[N] R;                // true recruit abundance (not density) by brood year
  // kelt survival
  vector[iter*N] epsilon_SS;           // process error in kelt survival by outmigration year
  vector<lower=0,upper=1>[iter*N] s_SS; // true kelt survival by outmigration year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                 // true p_HOS in all years (can == 0)
  matrix<lower=0>[N,N_age+iter] S_W_a; // true wild spawner abundance by age (if iter, max is plus-group)
  vector<lower=0>[N] S_W;              // true total wild spawner abundance
  vector[N] S_H;                       // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  vector[N_age-1] mu_alr_p[N_pop];     // population mean log ratio age distributions
  simplex[N_age] p[N];                 // true cohort (maiden) age distributions
  matrix<lower=0,upper=1>[N,iter ? N_age*2 : N_age] q; // true spawner or [maiden | kelt] age distns

  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;

  // S-R parameters including covariate effects
  alpha_Xbeta = alpha[pop] .* exp(rows_dot_product(X_alpha, beta_alpha[pop,]));
  Rmax_Xbeta = Rmax[pop] .* exp(rows_dot_product(X_Rmax, beta_Rmax[pop,]));

  // Log-ratio transform of pop-specific mean cohort age distributions
  for(j in 1:N_pop)
    mu_alr_p[j] = log(head(mu_p[j], N_age-1)) - log(mu_p[j,N_age]);

  // Calculate true total wild and hatchery spawners and spawner age distribution
  // and predict recruitment from brood year i
  for(i in 1:N)
  {
    int ii; // index into S_init and q_init
    // number of orphan (maiden) age classes <lower=0,upper=N_age>
    int N_orphan_age = max(N_age - max(pop_year[i] - min_age, 0), N_age);
    vector[N_orphan_age] q_orphan; // orphan (maiden) age distribution
    vector[N_age] alr_p; // alr(p[i,])

    // Spawners and age structure
    // use initial values for orphan age classes, otherwise use process model
    if(pop_year[i] <= max_age)
    {
      ii = (pop[i] - 1)*max_age + pop_year[i];
      q_orphan = append_row(sum(head(q_init[ii], N_age - N_orphan_age + 1)), 
                            tail(q_init[ii], N_orphan_age - 1));
    }
    
    // Multivariate Matt trick for within-pop, time-varying (maiden) age vectors
    alr_p = rep_vector(0,N_age);
    alr_p[1:(N_age-1)] = mu_alr_p[pop[i]] + sigma_p[pop[i],]' .* (L_p[pop[i]] * zeta_p[i,]');
    alr_p = exp(alr_p);
    p[i] = alr_p / sum(alr_p);

    // AR(1) recruitment process errors
    if(pop_year[i] == 1)
      epsilon_R[i] = zeta_R[i] * sigma_R[pop[i]] / sqrt(1 - rho_R[pop[i]]^2);
    else
      epsilon_R[i] = rho_R[pop[i]] * epsilon_R[i-1] + sigma_R[pop[i]] * zeta_R[i];

    if(iter)  // iteroparous
    {
      row_vector[N_age+1] S_M_a = rep_row_vector(0, N_age + 1); // true wild maiden spawners by age
      row_vector[N_age+1] S_K_a = rep_row_vector(0, N_age + 1); // true wild kelt spawners by age
      
      // AR(1) kelt survival process errors
      if(pop_year[i] == 1)
        epsilon_SS[i] = zeta_SS[i] * sigma_SS[pop[i]] / sqrt(1 - rho_SS[pop[i]]^2);
      else
        epsilon_SS[i] = rho_SS[pop[i]] * epsilon_SS[i-1] + sigma_SS[pop[i]] * zeta_SS[i];
      
      s_SS[i] = inv_logit(logit(mu_SS[pop[i]]) + dot_product(X_SS[i], beta_SS[pop[i],]) + epsilon_SS[i]);  

      // Maiden spawners
      for(a in 1:N_age)
      {
        if(pop_year[i] <= ages[a]) // use initial values
        {
          if(pop_year[i] == 1) // use [maiden | kelt] initial age dist
            S_M_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_iter_init[pop[i]][a];
          else // use maiden-only initial age dist
            S_M_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_orphan[a - (N_age - N_orphan_age)];
        }
        else // use recruitment process model
          S_M_a[a] = R[i-ages[a]]*p[i-ages[a],a]*(1 - age_F[a]*F_rate[i])*(1 - age_B[a]*B_rate_all[i]);
      }
      
      // Kelts
      if(pop_year[i] == 1) // use initial values
        S_K_a[2:] = S_init[ii]*(1 - p_HOS_all[i])*tail(q_iter_init[pop[i]], N_age)';
      else  // use recruitment process model (pool plus group and max maiden age)
      {
        S_K_a[2:] = append_col(head(S_W_a[i-1,], N_age - 1), sum(tail(S_W_a[i-1,], 2))) * s_SS[i-1]; 
        S_K_a = S_K_a .* (1 - age_F' * F_rate[i]) .* (1 - age_B' * B_rate_all[i]);
      }
      
      S_W_a[i,] = S_M_a + S_K_a;    
      S_W[i] = sum(S_W_a[i,]);
      q[i,] = append_col(S_M_a[:N_age], S_K_a[2:])/S_W[i]; // [maiden | kelt] spawner age distribution
    }
    else  // semelparous
    {
      for(a in 1:N_age)
      {
        if(pop_year[i] <= ages[a]) // use initial values
          S_W_a[i,a] = S_init[ii]*(1 - p_HOS_all[i])*q_orphan[a - (N_age - N_orphan_age)];
        else  // use recruitment process model
          S_W_a[i,a] = R[i-ages[a]]*p[i-ages[a],a]*(1 - age_F[a]*F_rate[i])*(1 - age_B[a]*B_rate_all[i]);
      }
      
      S_W[i] = sum(S_W_a[i,]);
      q[i,] = S_W_a[i,]/S_W[i]; // total spawner age distribution
    }

    // Hatchery-origin and total spawners
    S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
    S[i] = S_W[i] + S_H[i];

    // Recruitment
    R_hat[i] = SR(SR_fun, alpha_Xbeta[i], Rmax_Xbeta[i], S[i], A[i]);
    R[i] = R_hat[i] * exp(dot_product(X_R[i], beta_R[pop[i],]) + epsilon_R[i]);
  }
}

model {
  vector[N_B] log_B_take; // log of true broodstock take when B_take_obs > 0

  // Priors

  // recruitment
  alpha ~ lognormal(prior_alpha[1], prior_alpha[2]);
  Rmax ~ lognormal(prior_Rmax[1], prior_Rmax[2]);
  to_vector(beta_R) ~ normal(0,5);
  rho_R ~ gnormal(0,0.85,20); // mildly regularize to ensure stationarity
  sigma_R ~ normal(0,5);
  zeta_R ~ std_normal();   // total recruits: R ~ lognormal(log(R_hat), sigma_R)

  // kelt survival
  mu_SS ~ beta(prior_mu_SS[1], prior_mu_SS[2]);
  to_vector(beta_SS) ~ normal(0,3);
  sigma_SS ~ normal(0,3);
  rho_SS ~ gnormal(0,0.85,20); // mildly regularize to ensure stationarity
  zeta_SS ~ normal(0,1);

  // (maiden) recruit age structure
  mu_p ~ dirichlet(prior_mu_p);
  to_vector(sigma_p) ~ normal(0,5);
  for(j in 1:N_pop)
    L_p[j] ~ lkj_corr_cholesky(1);
  // age probs logistic MVN: 
  // alr_p[i,] ~ MVN(mu_alr_p[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ std_normal();

  // removals
  log_B_take = log(S_W_a[which_B,]*age_B) + logit(B_rate); 
  // implies B_take[i] = (S_W_a[i] * ageB) * B_rate[i] / (1 - B_rate[i])
  B_take_obs ~ lognormal(log_B_take, 0.05); // penalty to force pred and obs broodstock take to match 

  // initial spawners and (maiden) spawner age distribution
  // (accounting for amalgamation of q_init to q_orphan)
  S_init ~ lognormal(mu_S_init, sigma_S_init);
  {
    matrix[N_age,max_age*N_pop] q_init_mat;
    for(j in 1:size(q_init)) q_init_mat[,j] = q_init[j];
    target += sum((mu_q_init - 1) .* log(q_init_mat)); // q_init[i] ~ Dir(mu_q_init[,i])
  }

  // spawner observation error
  tau ~ gnormal(prior_tau[1], prior_tau[2], prior_tau[3]);

  // Observation model
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau[pop[which_S_obs]]); // obs total spawners
  n_H_obs ~ binomial(n_HW_obs, p_HOS); // observed counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q)); // obs wild spawner age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  corr_matrix[N_age-1] R_p[N_pop]; // correlation matrices of within-pop cohort log-ratio age distns
  vector[N] LL_S_obs;              // pointwise log-likelihood of total spawners
  vector[N_H] LL_n_H_obs;          // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_age_obs;          // pointwise log-likelihood of wild spawner age age frequencies
  vector[N] LL;                    // total pointwise log-likelihood

  for(j in 1:N_pop)
    R_p[j] = multiply_lower_tri_self_transpose(L_p[j]);

  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau[pop[which_S_obs]]);
  LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1, iter ? N_age*2 : N_age);
  LL_n_H_obs = rep_vector(0,N_H);
  for(i in 1:N_H)
    LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  LL = LL_S_obs + LL_n_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
