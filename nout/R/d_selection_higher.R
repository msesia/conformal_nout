# For Lehmann's alternatives we know that the outlier distribution g is monotone

#' d_selection_higher
#'
#'@description  It performs closed testing method with (higher) WMW local tests using an exact shortcut
#' relying on the increasing monotonicity of the Lehmann's alternatives.
#'
#' @param S_X   calibration score vector.
#' @param S_Y  test score vector.
#' @param S  selection set in the index test set. If \code{NULL} the entire test set is selected.
#' @param local_test  it can be either "wmw" for Wilcoxon sum-rank test or "higher" for higher order Wilcoxon sum-rank test.
#' @param k  order of the generalized Wilcoxon rank sum test. Classic Wilcoxon test corresponds to \eqn{k=1}.
#' @param alpha  significance level.
#' @param pvalue_only  logical value. If \code{TRUE}, only the global test is performed.
#' @param n_perm  minimum test sample size needed to use the asymptotic distribution of the test statistic.
#' @param B  number of replications to compute critical values and global *p*-value. Default value is 10^3.
#' @param critical_values  if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic.
#' @param seed  seed to ensure reproducible results.
#'
#' @return A list:
#' \itemize{
#' \item \code{lower_bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using the chosen local test.
#' \item \code{S}: a vector which is the selection set. If \code{NULL}, the selection set is the entire test set.
#' \item \code{global.pvalue}: a number which is the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null.
#' \item \code{selection.pvalue}: a number which is the *p*-value for the selected null. By default it is set equal to 1.
#' }
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#' X = runif(10)
#' Y = replicate(10, rg2(rnull=runif))
#' res1 = d_selection_higher(X, Y, local_test="WMW", n_perm=0, B=100)
#' res2 = d_selection_higher(X, Y, local_test="higher", k=2, S = c(1:7), n_perm=0, B=100)
d_selection_higher = function(S_X, S_Y, S=NULL, local_test="wmw", k=NULL, alpha=0.1, pvalue_only=FALSE, n_perm=0, B=10^3, critical_values=NULL, seed=123) {

  local_test <- tolower(local_test)
  stopifnot(local_test %in% c("wmw", "higher"))

  if(local_test=="wmw") {
    k <- 1
  } else { stopifnot(k>1 & k%%1==0) }

  if(k==1) local_test="wmw"

  m = as.double(length(S_X))
  n = as.double(length(S_Y))
  N = as.double(n+m)
  s = ifelse(is.null(S), n, length(S))

  Z = c(S_X, base::sort(S_Y, decreasing = F))

  if(!pvalue_only) {

    if(local_test=="wmw") { # Use Mann-Whitney test statistic (ranks computed in the calibration set only) and Tian et al.(2023) shortcut

      # Compute individual statistics for each test point
      S_Z = c(S_X, S_Y)
      R = stat.MW(Z=S_Z, m=m)
      # Compute all critical values for (m,k) from k in {1,...,n}
      crit = sapply(1:n, function(h) as.double(stats::qnorm(alpha, mean=m*h/2, sd = sqrt(m*h*(m+h+1)/12), lower.tail = F)))

      # Compute lower bound for S
      S_arg <- if (is.null(S)) NULL else as.integer(S)
      res <- sumSome::sumStatsPar(g = R, S = S_arg, alpha = alpha, cvs = crit)

      ## Compute p-value for the global null
      if(is.null(S)) {
        R.S = R[1:n]
      } else {
        R.S = R[S]
      }
      T.global.S = sum(R.S)

      pval.global = stats::pnorm(q=T.global.S, mean=m*s/2, sd = sqrt(m*s*(m+s+1)/12), lower.tail = F)
      d_S = res$TD

    } else { # for higher order WMW tests use our shortcut
      ## Find d
      Z = c(S_X, base::sort(S_Y, decreasing = F))

      # ----- Performance refactor (mathematically identical to the original) -----
      #
      # Original:
      #   R = sapply(n:1, function(h) stat.Tk(Z=Z[1:(m+h)], m=m, k=k))
      #   crit = compute.critical.values(...)
      #   T_wc = sapply(n:1, function(h) sum(R[[h]]))
      #   d = sum(cumsum(rev(T_wc) > rev(crit)) == 1:n)
      #
      # The original computes T_wc[h] and crit[h] for ALL h=1..n up front, then
      # asks "starting from h=n descending, how many consecutive h's reject?".
      # In practice the chain breaks at the first non-rejection, so we can do
      # this lazily.
      #
      # Two optimizations:
      #
      # (1) ONE rank() call instead of n.
      #     With Y sorted ascending, the excluded test points in Z[(m+h+1):(m+n)]
      #     are larger than every element in Z[1:(m+h)]. Removing the largest
      #     elements does not change anyone else's rank. Therefore:
      #         rank(Z[1:(m+h)]) == rank(Z)[1:(m+h)] for all h.
      #
      # (2) Lazy closed-testing loop with early termination.
      #     Iterate h from n downward. Compute T_wc[h] and crit[h] on the fly.
      #     Break at the first non-rejection (same decision as the cumsum/==
      #     1:n formulation; no statistical change).
      #
      # The asymptotic critical value path goes through get_crit_h (see
      # critical_values.R) which uses the vectorized fast_asymp_crit_Tk
      # (see WMW_test.R) instead of the slow asymptotic.moments.Tk's
      # interpreted sapply.

      R_full   <- base::rank(Z)
      R_Y_full <- R_full[(m + 1):(m + n)]

      d <- 0
      T_global <- NA_real_

      for (h in n:1) {
        R_h <- R_Y_full[1:h]
        T_h <- sum(stat.Tk_from_R(R_h, N_h = m + h, k = k))

        if (h == s) T_global <- T_h

        crit_h <- get_crit_h(m = m, h = h, k = k, alpha = alpha,
                             n_perm = n_perm, B = B,
                             critical_values = critical_values,
                             seed = seed,
                             local_test = local_test)

        # NOTE: this should be strictly larger!
        if (T_h > crit_h) {
          d <- d + 1
        } else {
          break
        }
      }

      # If we broke before reaching h=s, compute T_global once.
      # Happens only when s < n AND we failed to reject before h=s.
      if (is.na(T_global)) {
        R_s <- R_Y_full[1:s]
        T_global <- sum(stat.Tk_from_R(R_s, N_h = m + s, k = k))
      }

      T.global <- T_global

      # Compute p-value for the global null
      pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=s, local_test="higher", k=k, n_perm=n_perm, B=B, seed=seed)

      # Compute p-value for the selected null
      # NOTE: this calculation is missing
      pval.selection = 1

      h = n - d

      ## Find d_S
      if (!is.null(S)) {

        Y_S <- sort(S_Y[S], decreasing = TRUE)
        s = length(Y_S)
        Y_notS <- sort(S_Y[-S], decreasing = FALSE)

        d_S = 0
        for (i in 1:s) {
          for (j in max(0, (h - s + i - 1)):0) {
            ZZ <- c(S_X, Y_S[i:s], Y_notS[0:j])

            # Use the helper to avoid the rank+apply work inside stat.Tk
            R_ZZ <- base::rank(ZZ)
            l <- length(ZZ) - m
            R_test <- R_ZZ[(m + 1):(m + l)]
            T_wc <- sum(stat.Tk_from_R(R_test, N_h = m + l, k = k))

            crit = get_crit_h(m = m, h = l, k = k, alpha = alpha,
                              n_perm = n_perm, B = B,
                              critical_values = critical_values,
                              seed = seed,
                              local_test = local_test)

            notrejected = T_wc < crit

            if (notrejected) {
              break
            }
          }

          if (notrejected) {
            break
          }
          d_S = d_S + 1

        }

      } else {
        d_S = d
      }

    }
  }  else { # if pvalue_only==TRUE (only when S=NULL)

    Y.S = sort(S_Y, decreasing = F)
    ZZ = c(S_X, Y.S)

    R = stat.Tk(Z=ZZ, m=m, k=k)
    T.global = sum(R)
    pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=n, local_test="higher", k=k, n_perm=n_perm, B=B, seed=seed)
    d_S=0
  }

  out = list("lower.bound" = d_S,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)

  return(out)

}
