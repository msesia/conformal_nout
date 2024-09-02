
#' d_selection_G
#'
#' @description  It performs closed testing method with Shiraishi local test using an exact shortcut
#' valid when the outlier density is monotone, either increasing or decreasing.
#' 
#' @param S_X :  calibration score vector
#' @param S_Y : test score vector
#' @param S : selection set in the index test set. If \code{NULL} the entire test set is selected
#' @param k : order of the generalized Wilcoxon rank sum test. Classic Wilcoxon sum-rank test corresponds to \eqn{k=1}
#' @param g.hat : it can be either a character ("analytical") or a function denoting the outlier density.
#' If g.hat=="analytical" the test statistics are computed analytically without Monte Carlo estimation.
#' If \code{NULL}, the outlier density is estimated from the data
#' @param monotonicity : character indicating if the outlier density function is monotone increasing or decreasing or neither. Default value is \code{NULL}
#' @param fit_method : character value indicating the method to approximate the outlier distribution when argument \code{g.hat} is \code{NULL}. 
#' It can be either "beta_mix" or "mixmodel"
#' @param prop.F  : proportion of inliers used to estimate the inlier distribution in the process of estimating the outlier density.
#' Default value is 0.5
#' @param alpha : significance level
#' @param pvalue_only : logical value. If TRUE, only the global test is performed
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
#' @param B_MC : number of replications to compute the Shiraishi test statistic
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
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#'
#' X = runif(10)
#' Y = replicate(10, rg2(rnull=runif))
#' res = d_selection_G(X, Y, S = c(1:7), g.hat = g2, monotonicity="increasing", B=100)
d_selection_G <- function(S_X, S_Y, S=NULL, k=NULL, g.hat=NULL, monotonicity=NULL, fit_method="beta_mix", prop.F=0.5, alpha=0.1, pvalue_only=FALSE, n_perm=10, B=10^3, B_MC=10^3, seed=123){
  
  if(!is.null(monotonicity))
    stopifnot("Error: monotonicity must be either increasing, decreasing"= monotonicity%in%c("decreasing", "increasing"))
  
  n = as.double(length(S_Y))
  m = as.double(length(S_X))
  N = as.double(m+n)
  s = ifelse(is.null(S), n, length(S))
  S_Z = c(S_X, S_Y)
  
  # If the outlier distribution is unknown it is estimated from the data
  if(is.null(g.hat)){
    monotone = ifelse(is.null(monotonicity), FALSE, TRUE)
    m1 = round(prop.F*m)
    S_X1 = sample(S_X,m1)
    S_pooled = c(setdiff(S_X,X1), Y)
    g.hat = estimate_g(S_X1, S_pooled, method=fit_method, monotone=monotone)$pdf
  }
  
  if(is.null(monotonicity))
    res = d_G_cons(S_X=S_X, S_Y=S_Y, S=S, g.hat=g.hat, k=k, alpha=alpha, pvalue_only=pvalue_only, n_perm=n_perm, B=B, B_MC=B_MC, seed=seed)
  else
    res = d_G_monotone(S_X=S_X, S_Y=S_Y, S=S, g.hat=g.hat, k=k, alpha=alpha, pvalue_only=pvalue_only, n_perm=n_perm, B=B, B_MC=B_MC, seed=seed)

  return(res)
  
}







#' d_G_monotone
#'
#' @description  It performs closed testing method with Shiraishi local test using an exact shortcut
#' valid when the outlier density is monotone, either increasing or decreasing.
#' 
#' @param S_X :  calibration score vector
#' @param S_Y : test score vector
#' @param S : selection set in the index test set. If \code{NULL} the entire test set is selected
#' @param g.hat : it can be either a character ("analytical") or a function denoting the outlier density.
#' If g.hat=="analytical" the test statistics are computed analytically without Monte Carlo estimation.
#' @param decr : logical value indicating whether the outlier distribution is decreasing (TRUE)
#' or increasing (FALSE)
#' @param k : order of the LMPI test statistic to be specified when g.hat is "analytical"
#' @param alpha : significance level
#' @param pvalue_only : logical value. If TRUE, only the global test is performed
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic when
#' local_test is either "higher" or "fisher"
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
#' @param B_MC : number of replications to compute the Shiraishi test statistic
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
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#' m = 10; n=10;
#' X = runif(m)
#' Y = replicate(n, rg2(rnull=runif))
#' res = d_G_monotone(X, Y, S=c(1:7), g.hat=g2, decr=F, B=100)
d_G_monotone = function(S_X, S_Y, S=NULL, g.hat, decr=F, k=NULL, alpha=0.1, pvalue_only=FALSE, n_perm=0, B=10^3, B_MC = 10^3, seed=123) {
  
  m = length(S_X)
  n = length(S_Y)
  N = m+n
  
  if(decr){
    X = -S_X
    Y = -S_Y
  } else {
    X = S_X
    Y = S_Y
  }
  Z = c(X, base::sort(Y, decreasing = FALSE))
  R = base::rank(Z)
  d = 0
  
  
  if(!pvalue_only){
    
    # Find d
    for (i in n:1) {
      
      # Monte Carlo simulation of elementary test statistics
      if(is.character(g.hat)){
        if(g.hat=="analytical"){
          a_i = sapply(1:(i+m), function(h) (k+1)*k_mom_beta(a=h, b=m+i-h+1, k=k))
        } else{
          cat("Error: g.hat must be either a density function or the string analytical.")
        }
        
      } else {
        a_i = apply(replicate(B_MC, g.hat(sort(stats::runif(m+i)))) , 1, mean)
      }
      
      if(decr){
        R_wc = rank(-Z[1:(m+i)])[(m+1):(m+i)]
      } else {
        R_wc = R[(m+1):(m+i)]
      }
      
      stat_i = sum(a_i[R_wc])
      mu_i = i * mean(a_i)
      var_i = i * m * sum((a_i - mean(a_i))^2) / ((m + i) * ((m + i) - 1))
      d = d + (stat_i > mu_i + qnorm(1 - alpha) * sqrt(var_i))
      
      # compute global p-value
      if(i==n){
        pval.global = compute.global.pvalue(T.obs=stat_i, m=m, n=n, local_test="g", stats_G_vector=a_i,
                                            n_perm=n_perm, B=B, seed=seed)
      }
      
      
      if (d < n - i + 1) {
        break
      }
    }
    h = n - d
    
    # Find d_S
    if (!is.null(S)) {
      
      if(decr){ # if outlier density is decreasing change the sign of the scores
        S_Y.S <- sort(S_Y[S], decreasing = FALSE)
        s = length(S_Y.S)
        S_Y.notS <- sort(S_Y[-S], decreasing = TRUE)
        X = -S_X
        Y_S = -S_Y.S
        Y_notS = -S_Y.notS
        
      } else { # if outlier density is increasing don't change the sign of the scores
        
        S_Y.S <- sort(S_Y[S], decreasing = TRUE)
        s = length(S_Y.S)
        S_Y.notS <- sort(S_Y[-S], decreasing = FALSE)
        X = S_X
        Y_S = S_Y.S
        Y_notS = S_Y.notS
      }
      
      a_N = sapply((m + 1):(m + n), function(N) {
        apply(replicate(B, g.hat(sort(runif(N)))), 1, mean)
      })
      
      d_S = 0
      
      for (i in 1:s) {
        for (j in max(0, (h - s + i - 1)):0) {
          
          ZZ <- c(X, Y_S[i:s], Y_notS[0:j])
          if(decr){
            RR <- base::rank(-ZZ)
          } else{
            RR <- base::rank(ZZ)
          }
          
          
          k = length(ZZ) - m
          stat_k <- sum(a_N[[k]][RR[(m + 1):(m + k)]])
          mu_k = k * mean(a_N[[k]])
          var_k = k * m * sum((a_N[[k]] - mean(a_N[[k]]))^2) / ((m + k) * ((m + k) - 1))
          notrejected <- stat_k < mu_k + qnorm(1 - alpha) * sqrt(var_k)
          
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
  } else { # if pvalue_only==TRUE (only when S=NULL)
    
    if(is.character(g.hat)){
      if(g.hat=="analytical"){
        a_i = sapply(1:(n+m), function(h) (k+1)*k_mom_beta(a=h, b=m+n-h+1, k=k))
      } else{
        cat("Error: g.hat must be either a density function or the string analytical.")
      }
      
    } else {
      a_i = apply(replicate(B_MC, g.hat(sort(stats::runif(m+n)))) , 1, mean)
    }
    
    if(decr){
      R_wc = rank(-Z[1:(m+n)])[(m+1):(m+n)]
    } else {
      R_wc = R[(m+1):(m+n)]
    }
    
    stat_i = sum(a_i[R_wc])
    mu_i = n * mean(a_i)
    var_i = n * m * sum((a_i - mean(a_i))^2) / ((m + n) * ((m + n) - 1))
    
    # compute global p-value
    pval.global = compute.global.pvalue(T.obs=stat_i, m=m, n=n, local_test="g", stats_G_vector=a_i,
                                        n_perm=n_perm, B=B, seed=seed)
    d_S=0
  }
  out = list("lower.bound" = d_S,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)
  
  return(out)
  
}









#' d_G_cons
#' @description  It performs closed testing method with Shiraishi local test using a shortcut
#' approximating the lower bound for the number of outliers when the outlier density is not monotone
#' 
#' @param S_X :  calibration score vector
#' @param S_Y : test score vector
#' @param S : selection set in the index test set. If \code{NULL} the entire test set is selected
#' @param g.hat : it can be either a character ("analytical") or a function denoting the outlier density.
#' If g.hat=="analytical" the test statistics are computed analytically withuout Monte Carlo estimation.
#' @param k : order of the LMPI test statistic to be specified when g.hat is "analytical"
#' @param alpha : significance level
#' @param pvalue_only : logical value. If TRUE, only the global test is performed
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
#' @param B_MC : number of replications to compute the Shiraishi test statistic
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
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#' m = 10; n=10;
#' X = runif(m)
#' Y = replicate(n, rg2(rnull=runif))
#' res = d_G_cons(X, Y, g.hat=g2, k=1, B=100)
d_G_cons = function(S_X, S_Y, S=NULL, g.hat, k=NULL, alpha=0.1, pvalue_only=FALSE, n_perm=10, B=10^3, B_MC=10^3, seed=123){

  m = as.double(length(S_X))
  n = as.double(length(S_Y))
  N = as.double(n+m)
  s = ifelse(is.null(S), n, length(S))

  if(!pvalue_only){
    tentative.d = 0; l = n
    cont = TRUE
    pval.global = 1

    while(cont == TRUE & l>0){
      if(is.character(g.hat)){
        if(g.hat=="analytical"){
          stats_G = sapply(1:(l+m), function(h) (k+1)*k_mom_beta(a=h, b=m+l-h+1, k=k))
        } else{
          cat("Error: g.hat must be either a density function or the string analytical.")
        }

      } else {
        stats_G = apply(replicate(B_MC, sapply(X=sort(stats::runif(m+l)), FUN=g.hat)), 1, mean)
      }

      if(is.null(S)){

        Rx = sapply(1:n, function(i) rank(c(S_X, S_Y[i]))[m+1]-1)
        range_test_ranks_n = (min(Rx)+1):(max(Rx)+l)
        R = sort(stats_G[range_test_ranks_n], decreasing=F)[1:l]

      } else {
        Rx = sapply(1:n, function(i) rank(c(S_X, S_Y[i]))[m+1]-1)
        Rx.S = Rx[S]
        if(l>=(n-s+1)){
          range_test_ranks_supS = (min(Rx)+1):(max(Rx)+l)
          range_test_ranks_S = range_test_ranks_supS
        } else {
          range_test_ranks_subS = (min(Rx.S)+1):(max(Rx.S)+l)
          range_test_ranks_S = range_test_ranks_subS
        }
        # R = stats_G[range_test_ranks_S]
        R = sort(stats_G[range_test_ranks_n], decreasing = F)[1:l]
      }

      T_wc = sum(R)
      crit = asymptotic.critical.G(m=m, n=l, stats_G_vector=stats_G, alpha=alpha)

      if(l==s){
        T.global = T_wc
        pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=s, local_test="g", stats_G_vector=stats_G,
                                            n_perm=n_perm, B=B, seed=seed)
      } else {
        pval.global = pval.global
      }

      if(T_wc > crit){ # NOTE: this should be strictly larger!
        tentative.d = tentative.d+1
        l=l-1
      } else{
        tentative.d = tentative.d
        cont=FALSE
      }
    }

    d = ifelse(tentative.d-n+s>0, tentative.d-n+s, 0)

  } else {
    Rx = sapply(1:n, function(i) rank(c(S_X, S_Y[i]))[m+1]-1)
    range_test_ranks_n = (min(Rx)+1):(max(Rx)+n)
    R = sort(stats_G[range_test_ranks_n], decreasing=F)[1:n]

    T.global = sum(R)
    pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=n, local_test="g", stats_G_vector=stats_G,
                                        n_perm=n_perm, B=B, seed=seed)
    d=0
  }

  out = list("lower.bound" = d,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)

}
