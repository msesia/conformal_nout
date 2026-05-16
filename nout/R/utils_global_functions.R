


#' perm.crit.T
#'
#' @description It computes the permuted \eqn{(1-\alpha)}-quantile of a chosen test statistic.
#'
#' @param m  calibration sample size.
#' @param n  test sample size.
#' @param local_test  local test to be used in the closed testing procedure. Default value is Wilcoxon sum-rank test.
#' @param k  order of the higher order Wilcoxon sum-rank test statistic.
#' @param stats_G_vector  vector of elementary test statistics to perform the test in Shiraishi (1985).
#' @param alpha  significance level. Default value is set equal to 0.1.
#' @param B  number of permutations.
#' @param seed  seed to ensure reproducible results.
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of the chosen test statistic obtained via permutation.
#'
#'
perm.crit.T <- function(m, n, local_test="wmw", k=NULL, stats_G_vector=NULL, alpha=0.1, B=10^3, seed=123){
  set.seed(seed)

  local_test = tolower(local_test)
  stopifnot(local_test %in% c("wmw", "higher", "fisher", "g"))

  if(local_test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local_test = "wmw"
  }
  if(local_test=="wmw") k=1

  T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    Z = stats::runif(N)
    if(local_test=="fisher"){
      T = sum(stat.Fisher(Z=Z, m=m))
    } else if(local_test=="g"){
      T = sum(stat.G(Z=Z, m=m, stats_G_vector=stats_G_vector))
    } else if(local_test=="higher" || local_test=="wmw"){
      T = sum(stat.Tk(Z=Z, k=k, m=m))
    }
  }

  # Empirical quantile with finite-sample adjustment
  idx = ceiling(B*(1-alpha)*(1+1/B))
  if(idx <= B) {
    crit = sort(as.vector(T.v))[idx]
  } else {
    crit = Inf
  }

  return(crit)
}



#' compute.perm.pval
#'
#' @description It computes the *p*-value for the global null obtained via permutation using a chosen test statistic.
#'
#' @param T.obs  observed value of the chosen test statistic.
#' @param local_test  local test to be used in the closed testing procedure. Default value is Wilcoxon sum-rank test.
#' @param m  calibration sample size.
#' @param n  test sample size.
#' @param k  order of the higher order Wilcoxon sum-rank test statistic.
#' @param stats_G_vector  vector of elementary test statistics to perform the test in Shiraishi (1985).
#' @param B  number of permutations.
#' @param seed  seed to ensure reproducible results.
#'
#'
#' @return It returns the *p*-value for the global null obtained via permutation using the chosen test statistic.
#'
#'
compute.perm.pval <- function(T.obs, local_test="wmw", m, n, k=NULL, stats_G_vector, B=10^3, seed=123) {

  set.seed(seed)

  local_test = tolower(local_test)
  stopifnot(local_test %in% c("wmw", "higher", "fisher", "g"))

  if(local_test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local_test = "wmw"
  }
  if(local_test=="wmw") k=1

  T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    Z = stats::runif(N)
    if(local_test=="fisher")
      T = sum(stat.Fisher(Z=Z, m=m))
    else if(local_test=="g")
      T = sum(stat.G(Z=Z, m=m, stats_G_vector=stats_G_vector))
    else if(local_test=="higher" || local_test=="wmw")
      T = sum(stat.Tk(Z=Z, k=k, m=m))
    return(T)
  }

  # Compute the permutation p-value
  pval = (1+sum(T.v >= T.obs)) / (1 + length(T.v))

  return(pval)
}



#' compute.global.pvalue
#'
#' @description It computes the *p*-value for the global null according to the chosen
#' test statistic. The *p*-value is computed via permutation if either the calibration
#' sample size or the test sample size is smaller than \code{n_perm}.
#' Otherwise, it is computed using the asymptotic distribution.
#'
#' @param T.obs  observed value of the chosen test statistic.
#' @param m  calibration sample size.
#' @param n  test sample size.
#' @param local_test  local test to be used in the closed testing procedure. Default value is Wilcoxon sum-rank test.
#' @param k  order of the higher order Wilcoxon sum-rank  test statistic.
#' @param stats_G_vector  vector of elementary test statistics to perform the test in Shiraishi (1985).
#' @param n_perm  if \eqn{min(m,n)\leq n_perm} the *p*-value for the global null will be computed via permutation. Default value is 0.
#' @param B  number of permutations.
#' @param seed  seed to ensure reproducible results.
#'
#'
#' @return A number, the *p*-value for the global null according to the chosen test statistic.
#'
#'
compute.global.pvalue <- function(T.obs, m, n, local_test="wmw", k=NULL, stats_G_vector, n_perm=0, B=100, seed=321) {

  local_test = tolower(local_test)
  stopifnot(local_test %in% c("wmw", "higher", "fisher", "g"))

  if(local_test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local_test = "wmw"
  }
  if(local_test=="wmw") k=1

  # permutation p-value for the global null if the sample size is small
  if(min(m,n)<=n_perm){
    if(local_test=="fisher")
      pval.perm = compute.perm.pval(T.obs=T.obs, local_test="fisher", m=m, n=n, k=k, B=B, seed=seed)
    else if(local_test=="g")
      pval.perm = compute.perm.pval(T.obs=T.obs, local_test="g", m=m, n=n, k=k, stats_G_vector=stats_G_vector, B=B, seed=seed)
    else if(local_test=="higher" || local_test=="wmw")
      pval.perm = compute.perm.pval(T.obs=T.obs, local_test="higher", m=m, n=n, k=k, B=B, seed=seed)
  }
  # Otherwise, use the asymptotic approximation to compute an approximate p-value for the global null
  else {
    if(local_test=="fisher")
      pval.perm = asymptotic.pvalue.Fisher(m=m, n=n, T.obs=T.obs)
    else if(local_test=="g")
      pval.perm = asymptotic.pvalue.G(m=m, n=n, T.obs=T.obs, stats_G_vector = stats_G_vector)
    else if(local_test=="higher" || local_test=="wmw")
      pval.perm = asymptotic.pvalue.Tk(m=m, n=n, k=k, T.obs=T.obs)
  }
  return(pval.perm)
}



#' compute.critical.values
#'
#' @description It computes the vector of critical values for a chosen test statistic
#' at significance level \eqn{\alpha}, one for each integer between 1 and the test size.
#'
#' @param m  calibration size.
#' @param n  test size.
#' @param local_test  local test to be used in the closed testing procedure. Default value is Wilcoxon sum-rank test.
#' @param alpha  significance level.
#' @param k  order of the higher order Wilcoxon sum-rank test statistic.
#' @param stats_G_vector  vector of elementary test statistics to perform the test in Shiraishi (1985).
#' @param n_perm  if \eqn{min(m,n)\leq n_perm} critical values will be computed via permutation. Default value is 0.
#' @param B  number of permutation to compute critical values. Default value is 10^3.
#' @param critical_values  if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic.
#' @param seed  seed to ensure reproducible results.
#'
#'
#' @return A vector of numeric values that are the critical values at significance
#' level \eqn{\alpha} for a test statistic chosen among higher order Wilcoxon sum-rank,
#' Fisher or Shiraishi statistics.
#'
compute.critical.values <- function(m, n, local_test="wmw", alpha, k=NULL, stats_G_vector, n_perm=0, B=10^3, critical_values=NULL, seed=123){

  local_test = tolower(local_test)
  stopifnot(local_test %in% c("wmw", "higher", "fisher", "g"))

  if(local_test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local_test = "wmw"
  }
  if(local_test=="wmw") k=1

  crit = sapply(1:n, function(h) {
    # For small values of m and n compute critical values via permutation
    if(min(m,h)<=n_perm) {
      found.value = FALSE
      # In order to avoid repeating computation of precomputed critical values that are saved in "tables" folder
      if(!is.null(critical_values)) {
        if(length(critical_values)>=h) {
          critical.value = critical_values[h]
          found.value = TRUE
        }
      }
      if(!found.value) {
        #cat(sprintf("Running permutations...\n"))

        critical.value = as.double(perm.crit.T(m=m, n=h, k=k, local_test=local_test, alpha=alpha, B=B, seed=seed))
      }
    }
    # For large values of m or n compute critical values using the asymptotic approximation of the test statistic
    else {
      if(local_test=="fisher"){
        critical.value = as.double(asymptotic.critical.Fisher(m=m, n=h, alpha=alpha))
      } else if(local_test=="g"){
        critical.value = as.double(asymptotic.critical.G(m=m, n=h, stats_G_vector=stats_G_vector, alpha=alpha))
      } else if(local_test=="higher" || local_test=="wmw"){
        critical.value = as.double(asymptotic.critical.Tk(m=m, n=h, k=k, alpha=alpha))
      }
    }
    return(as.double(critical.value))
  })
  return(crit)
}




#' compute.1critical.value
#'
#' @description It computes the critical value for a chosen test statistic
#' at significance level \eqn{\alpha}.
#'
#' @param m  calibration size.
#' @param n  test size.
#' @param local_test  local test to be used in the closed testing procedure. Default value is Wilcoxon sum-rank test.
#' @param alpha  significance level.
#' @param k  order of the higher order Wilcoxon sum-rank test statistic.
#' @param stats_G_vector  vector of elementary test statistics to perform the test in Shiraishi (1985).
#' @param n_perm  if \eqn{min(m,n)\leq n_perm} critical values will be computed via permutation. Default value is 0.
#' @param B  number of permutation to compute critical values. Default value is 10^3.
#' @param critical_values  if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic.
#' @param seed  seed to ensure reproducible results.
#'
#'
#' @return A numeric value that is the critical value at significance level \eqn{\alpha} for a test statistic chosen among
#' higher order Wilcoxon sum-rank, Fisher or Shiraishi statistics.
#'
#'
compute.1critical.value <- function(m, n, local_test="wmw", alpha, k=NULL, stats_G_vector, n_perm=0, B=10^3, critical_values=NULL, seed=123){

  local_test = tolower(local_test)
  stopifnot(local_test %in% c("wmw", "higher", "fisher", "g"))

  if(local_test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local_test = "wmw"
  }
  if(local_test=="wmw") k=1

  # For small values of m and n compute critical values via permutation
  if(min(m,n)<=n_perm) {
    found.value = FALSE
    # In order to avoid repeating computation of precomputed critical values that are saved in "tables" folder
    if(!is.null(critical_values)) {
      if(length(critical_values)>=n) {
        critical.value = critical_values[n]
        found.value = TRUE
      }
    }
    if(!found.value) {
      #cat(sprintf("Running permutations...\n"))

      critical.value = as.double(perm.crit.T(m=m, n=n, k=k, local_test=local_test, alpha=alpha, B=B, seed=seed))
    }
  } else { # For large values of m or n compute critical values using the asymptotic approximation of the test statistic
    if(local_test=="fisher"){
      critical.value = as.double(asymptotic.critical.Fisher(m=m, n=n, alpha=alpha))
    } else if(local_test=="g"){
      critical.value = as.double(asymptotic.critical.G(m=m, n=n, stats_G_vector=stats_G_vector, alpha=alpha))
    } else if(local_test=="higher" || local_test=="wmw"){
      critical.value = as.double(asymptotic.critical.Tk(m=m, n=n, k=k, alpha=alpha))
    }
  }

  return(critical.value)
}

#' get_crit_h
#'
#' @description Critical value for one test sample size \eqn{h}, with the same
#' dispatch as \code{compute.1critical.value} but using the vectorized
#' \code{fast_asymp_crit_Tk} for the asymptotic branch. Used by
#' \code{d_selection_higher}'s closed-testing loop to avoid the slow
#' \code{asymptotic.moments.Tk} sapply on each iteration.
#'
#' Dispatch (matches \code{compute.1critical.value} exactly):
#' \itemize{
#'   \item If \eqn{\min(m, h) \le n_{perm}}: delegate to
#'         \code{compute.1critical.value}, which itself reads from
#'         \code{critical_values[h]} if available or runs permutation
#'         otherwise. Behavior is preserved EXACTLY (including seeding and
#'         foreach evaluation), so the permutation branch is unchanged.
#'   \item Otherwise: use the vectorized \code{fast_asymp_crit_Tk}, which is
#'         mathematically identical to \code{asymptotic.critical.Tk}.
#' }
#'
#' This function currently only optimizes the \code{"higher"} (and equivalently
#' \code{"wmw"}) branch. The \code{"fisher"} and \code{"g"} branches still
#' delegate fully to \code{compute.1critical.value}.
#'
#' @param m                calibration size.
#' @param h                test sample size at which to evaluate.
#' @param k                order of the higher-order WMW statistic. Required
#'                         for \code{local_test = "higher"}.
#' @param alpha            significance level.
#' @param n_perm           threshold below which the permutation branch is used.
#' @param B                number of permutations (passed through).
#' @param critical_values  optional precomputed table.
#' @param seed             RNG seed (passed through).
#' @param local_test       "wmw", "higher", "fisher", or "g". Default "higher".
#'
#' @return A numeric scalar: the critical value at \eqn{h} for level \eqn{\alpha}.
#'
#' @keywords internal
get_crit_h <- function(m, h, k, alpha, n_perm, B, critical_values, seed,
                       local_test = "higher") {
  if (min(m, h) <= n_perm) {
    return(as.double(compute.1critical.value(
      m = m, n = h, local_test = local_test,
      alpha = alpha, k = k, n_perm = n_perm, B = B,
      critical_values = critical_values, seed = seed
    )))
  }
  ## Asymptotic branch. Only the higher/wmw case has a vectorized fast path;
  ## others delegate to the package.
  if (local_test %in% c("higher", "wmw")) {
    return(fast_asymp_crit_Tk(m = m, h = h, k = k, alpha = alpha))
  }
  as.double(compute.1critical.value(
    m = m, n = h, local_test = local_test,
    alpha = alpha, k = k, n_perm = n_perm, B = B,
    critical_values = critical_values, seed = seed
  ))
}

