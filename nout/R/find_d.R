#' find_d
#'
#' @param X : calibration score vector
#' @param Y : test score vector
#' @param local_test : local test to be used in the closed testing procedure.
#' It can be either "wmw" for Wilcoxon sum-rank test, "higher" for higher order Wilcoxon sum-rank tests,
#' "fisher" for Fisher's combination test, "g" for the test by Shiraishi (1985), "simes" for Simes' test or
#' "storey" for Simes' test using Storey's estimator for the proportion of true null hypotheses.
#' @param S : selection set in the index test set
#' @param k : positive integer indicating the order of generalized Wilcoxon sum-rank test. Default value is \code{NULL}
#' @param monotonicity : character indicating if the outlier density function is monotone increasing or decreasing or neither. Default value is \code{NULL}
#' @param g.hat : it denotes the outlier density. If \code{NULL}, it is estimated from the data
#' @param fit_method : character value indicating the method to approximate the outlier distribution when argument \code{g.hat} is \code{NULL}. 
#' It can be either "beta_mix" or "mixmodel"
#' @param prop.F  : proportion of inliers used to estimate the inlier distribution in the process of estimating the outlier density.
#' Default value is 0.5
#' @param pvalue_only : logical value. If TRUE, only the global test is performed
#' @param alpha : significance level
#' @param lambda : parameter to be specified when computing Storey's estimator. Default value is 0.5
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
#' @param B_MC : number of replications to compute the Shiraishi test statistic
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param seed : seed to ensure reproducible results
#' 
#' @return A list:
#' \itemize{
#' \item \code{lower_bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using the chosen local test
#' \item \code{S}: a vector which is the selection set. If \code{NULL}, the selection set is the entire test set
#' \item \code{global.pvalue}: a number which is the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null
#' \item \code{selection.pvalue}: a number which is the *p*-value for the selected null. By default it is set equal to 1
#' }
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#'
#' X = runif(10)
#' Y = replicate(10, rg2(rnull=runif))
#' res = find_d(X, Y, local_test="higher", k=3, B=100)
#' res = find_d(X, Y, local_test="g", g.hat = g2, monotonicity="increasing", B=100)
#' 
find_d = function(X, Y, local_test = "wmw", S=NULL, k=NULL, monotonicity=NULL, g.hat=NULL, fit_method="beta_mix", prop.F=0.5, pvalue_only=FALSE, alpha=0.1, lambda=0.5, n_perm=0, B=10^3, B_MC=10^3, critical_values=NULL, seed=123){

  local_test = tolower(local_test)
  stopifnot(local_test %in% c("wmw", "higher", "fisher", "g", "simes", "storey"))

  if(local_test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local_test = "wmw"
  }

  if(local_test=="wmw") k=1


  if(local_test=="wmw" || local_test=="higher"){

    res = d_selection_higher(X, Y, S=S,local_test=local_test, k=k, alpha=alpha, pvalue_only=pvalue_only, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed )

  } else if(local_test=="fisher"){

    res = d_selection_fisher(X, Y, S=S, alpha=alpha, pvalue_only=pvalue_only, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed )

  } else if(local_test=="g"){

    res = d_selection_G(X, Y, S=S, k=k, g.hat=g.hat, monotonicity=monotonicity, fit_method=fit_method, prop.cal=prop.cal, alpha=alpha, pvalue_only=pvalue_only, n_perm=n_perm, B=B, B_MC=B_MC, seed=seed)

  } else if(local_test=="simes"){

    res = d_selection_simes(X, Y, S=S, alpha=alpha)

  } else if(local_test=="storey"){

    res = d_selection_storey(X, Y, S=S, alpha=alpha, lambda=lambda)

  }

  return(res)

}



