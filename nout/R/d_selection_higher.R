# For Lehmann's alternatives we know that the outlier distribution g is monotone

#' d_selection_higher
#'
#'@description  It performs closed testing method with (higher) WMW local tests using an exact shortcut
#' relying on the increasing monotonicity of the Lehmann's alternatives.
#' 
#' @param S_X :  calibration score vector
#' @param S_Y : test score vector
#' @param S : selection set in the index test set
#' @param local_test : it can be either "wmw" for Wilcoxon rank sum test or "higher" for higher order Wilcoxon rank sum test
#' @param k : order of the generalized Wilcoxon rank sum test. Classic Wilcoxon test corresponds to \eqn{k=1}
#' @param alpha : significance level
#' @param pvalue_only : logical value. If TRUE, only the global test is performed
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic when
#' local_test is either "higher" or "fisher"
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param seed : seed to ensure reproducible results
#'
#' @return A list:
#' \itemize{
#' \item \code{lower.bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using the chosen local test
#' \item \code{S}: the selection set, i.e., the selected subset of the test indices
#' \item \code{global.pvalue}: the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null
#' \item \code{selection.pvalue}: *p*-value for the selected null
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
d_selection_higher = function(S_X, S_Y, S=NULL, local_test="wmw", k=NULL, alpha=0.1, pvalue_only=FALSE, n_perm=0, B=10^3, critical_values=NULL, seed=123){
  
  local_test=tolower(local_test)
  stopifnot(local_test %in% c("wmw", "higher"))
  
  if(local_test=="wmw") {
    k=1
  } else { stopifnot(k>1 & k%%1==0) }
  
  if(k==1) local_test="wmw"
  
  m = as.double(length(S_X))
  n = as.double(length(S_Y))
  N = as.double(n+m)
  s = ifelse(is.null(S), n, length(S))
  
  Z = c(S_X, base::sort(S_Y, decreasing = F))
  
  if(!pvalue_only){
    
    if(local_test=="wmw"){ # Use Mann-Whitney test statistic (ranks computed in the calibration set only) and Tian et al.(2023) shortcut
      
      # Compute individual statistics for each test point
      S_Z = c(S_X, S_Y)
      R = stat.MW(Z=S_Z, m=m)
      
      # Compute all critical values for (m,k) from k in {1,...,n}
      crit =sapply(1:n, function(h) as.double(stats::qnorm(alpha, mean=m*h/2, sd = sqrt(m*h*(m+h+1)/12), lower.tail = F)))
      
      # Compute lower bound for S
      res = sumSome::sumStatsPar(g = R, S = S, alpha = alpha, cvs = crit)
      
      ## Compute p-value for the global null
      if(is.null(S)){
        R.S = R[1:n] 
      } else {
        R.S = R[S]
      }
      T.global.S = sum(R.S)
      
      pval.global = stats::pnorm(q=T.global.S, mean=m*s/2, sd = sqrt(m*s*(m+s+1)/12), lower.tail = F)
      d_S = res$TD
      
      if(is.null(S)){
        S=NULL
      }
      
    } else { # for higher order WMW tests use our shortcut
      ## Find d
      Z = c(X,base::sort(S_Y, decreasing = F))
      
      # Closed-testing shortcut: sort the test points based on their individual statistics
      # For each k in {1,...,n} consider the worst-case subset of test points with cardinality k
      R = sapply(n:1, function(h) stat.Tk(Z=Z[1:(m+h)], m=m, k=k))
      
      # Compute all critical values for (m,k) from k in {1,...,n}
      crit = as.double(compute.critical.values(m=m, n=n, local_test=local_test,
                                               alpha=alpha, k=k, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed))
      
      T_wc = sapply(n:1, function(h) sum(R[[h]]))
      
      # Compare the worst-case statistics to the critical values for k in {n,...,1}, starting from the max cardinality
      d = as.double(sum(cumsum(rev(T_wc) > rev(crit)) == 1:n)) # NOTE: this should be strictly larger!
      
      # Compute p-value for the global null
      T.global = T_wc[s]
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
            R = stat.Tk(Z=ZZ, m=m, k=k)
            l = length(ZZ) - m
            crit = as.double(compute.1critical.value(m=m, n=l, local_test=local_test,
                                                     alpha=alpha, k=k, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed))
            
            T_wc = sum(R)
            
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




