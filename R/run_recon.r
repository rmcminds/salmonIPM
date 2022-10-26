#' Perform run reconstruction on brood table data
#'
#' @param fish_data Data frame that includes the following columns, 
#' in no particular order except where noted:
#' * `pop`  Numeric, factor or character population ID.
#' * `year`  Numeric variable giving the year the fish spawned (i.e., the brood year).
#' * `A` Spawning habitat size (either stream length or area).
#' Will often be time-invariant within a population, but need not be.
#' * `S_obs`  Total number (not density) of all wild and hatchery-origin spawners.
#' * `n_age[min_age]_obs...n_age[max_age]_obs`  Multiple columns of
#'   observed spawner age frequencies (i.e., counts), where `[min_age]` and
#'   `[max_age]` are the numeral age in years (total, not ocean age) of the
#'   youngest and oldest spawners, respectively. Complex age structures such as
#'   Gilbert-Rich or maiden-kelt are not supported.
#' * `n_W_obs`  Observed frequency of natural-origin spawners.
#' * `n_H_obs`  Observed frequency of hatchery-origin spawners.
#' * `F_rate`  Total harvest rate (proportion) of natural-origin fish.
#' * `B_take_obs`  Number of adults taken for hatchery broodstock.
#' @inheritParams salmonIPM
#' 
#' @return A data frame with the following columns, some of which are simply 
#' replicated from `fish_data`:
#' 
#' * `pop`  See above.
#' * `year`  See above.
#' * `A`  See above.
#' * `S_obs`  See above.
#' * `q_age[minAge]_obs...q_age[maxAge]_obs`  Multiple columns of spawner age 
#' proportions corresponding to the frequencies in `fish_data`.
#' * `p_age_minAge...p_age_maxAge`  Multiple columns of recruit age 
#' proportions by brood year.
#' * `p_HOS_obs`  Proportion of hatchery-origin spawners.
#' * `F_rate`  See above.
#' * `B_take_obs`  See above.
#' * `R_obs`  Total natural-origin recruits from the brood year in each row.
#' 
#' @export

run_recon <- function(fish_data, F_ages = NULL, B_ages = NULL)
{
  # Basic objects and data dimensions
  fish_data <- as.data.frame(fish_data)
  N <- nrow(fish_data)
  for(i in names(fish_data)) assign(i, fish_data[[i]])
  n_age_obs <- fish_data[, grep("n_age", names(fish_data))]
  N_age <- ncol(n_age_obs)
  n_obs <- rowSums(n_age_obs)
  ages <- as.numeric(substring(names(n_age_obs), 6, 6))
  q_obs <- sweep(n_age_obs, 1, n_obs, "/")
  for(i in unique(pop))
    if(any(pop == i & n_obs == 0))
      q_obs[pop == i & n_obs == 0,] <- rep(colMeans(q_obs[pop == i & n_obs > 0,]),
                                           each = sum(pop == i & n_obs==0))
  names(q_obs) <- gsub("n", "q", names(q_obs)) 
  p_HOS_obs <- ifelse(n_H_obs + n_W_obs > 0, n_H_obs / (n_H_obs + n_W_obs), 0)
  
  # Fishery and broodstock age-selectivity
  if(is.null(F_ages)) F_ages <- ages
  if(length(F_ages) == N_age) {
    if(is.logical(F_ages)) F_ages <- ages[F_ages]
    if(is.integer(F_ages) & all(F_ages %in% 0:1)) F_ages <- ages[as.logical(F_ages)]
  }
  
  if(is.null(B_ages)) B_ages <- ages
  if(length(B_ages) == N_age) {
    if(is.logical(B_ages)) B_ages <- ages[B_ages]
    if(is.integer(B_ages) & all(B_ages %in% 0:1)) B_ages <- ages[as.logical(B_ages)]
  }
  which_B_age <- B_ages - min(ages) + 1
  
  B_rate <- B_take_obs / (S_obs * (1 - p_HOS_obs) * rowSums(q_obs[,which_B_age]) + B_take_obs)

  # Reconstruct recruits and age structure by brood year and return results
  R_a <- matrix(NA, N, length(ages))
  for(i in 1:N)
    for(j in 1:length(ages)) {
      a <- ages[j]
      if(year[i] + a <= max(year[pop==pop[i]])) {
        F_eff <- ifelse(a %in% F_ages, F_rate[i+a], 0)
        B_eff <- ifelse(a %in% B_ages, B_rate[i+a], 0)
        R_a[i,j] <- S_obs[i+a] * (1 - p_HOS_obs[i+a]) * q_obs[i+a,j] / ((1 - B_eff) * (1 - F_eff))
      }
    }

  R_obs <- rowSums(R_a)
  p_obs <- sweep(R_a, 1, R_obs, "/")
  colnames(p_obs) <- gsub("q", "p", names(q_obs))
  rr_dat <- cbind(pop = pop, A = A, year = year, S_obs = S_obs, R_obs = R_obs, q_obs, p_obs, 
                  p_HOS_obs = p_HOS_obs, F_rate = F_rate, B_take_obs = B_take_obs)
  return(rr_dat)
}
