#' d_selection_simes
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test.
#'
#' @param S_X : score vector of calibration observations.
#' @param S_Y : score vector of test observations.
#' @param S : selection set in the index test set. If \code{NULL} the entire test set is selected.
#' @param alpha : significance level of the local test. Default value is set equal to 0.1.
#' @param pvalue_only : logical value. If \code{TRUE}, only the global test is performed.
#'
#' @return A list:
#' \itemize{
#' \item \code{lower_bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using the chosen local test.
#' \item \code{S}: a vector which is the selection set. If \code{NULL}, the selection set is the entire test set.
#' \item \code{global.pvalue}: a number which is the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null.
#' \item \code{selection.pvalue}: a number which is the *p*-value for the selected null.
#' }
#'
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_selection_simes(Sx, Sy, S=3)
#' d_selection_simes(Sx, Sy, S=c(3, 7:13))
#' d_selection_simes(Sx, Sy)
#'
d_selection_simes = function(S_X, S_Y, S=NULL, alpha = 0.1, pvalue_only=FALSE){
  n = length(S_Y)
  m = length(S_X)
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  hom = hommel::hommel(pval)

  # Compute p-value for the global null
  pval.global = hommel::localtest(hom)

  if(!pvalue_only){
    if(is.null(S)){

      # Lower bound
      d = hommel::discoveries(hom, alpha = alpha)
      # Compute p-value for the selected null
      pval.selection = pval.global

    } else {

      # Lower bound
      d = hommel::discoveries(hom, ix=S, alpha = alpha)
      # Compute p-value for the selected null
      pval.selection = hommel::localtest(hom, ix=S)
    }
  } else {
    pval.selection = 1
    d = 0
  }



  out = list("lower.bound" = d,
             "global.p.value" = pval.global,
             "S"=S,
             "selection.p.value" = pval.selection)

  return(out)
}



#' d_selection_storey
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test with Storey estimator for the proportion of true null hypotheses.
#'
#' @param S_X : score vector of calibration observations.
#' @param S_Y : score vector of test observations.
#' @param S : selection set in the index test set. If \code{NULL} the entire test set is selected.
#' @param alpha : significance level of the local test. Default value is set equal to 0.1.
#' @param lambda : parameter involved in the computation of Storey estimator. Default value is set equal to 0.5.
#' @param pvalue_only : logical value. If \code{TRUE}, only the global test is performed.
#'
#'
#' @return A list:
#' \itemize{
#' \item \code{lower_bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using the chosen local test.
#' \item \code{S}: a vector which is the selection set. If \code{NULL}, the selection set is the entire test set.
#' \item \code{global.pvalue}: a number which is the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null.
#' \item \code{selection.pvalue}: a number which is the *p*-value for the selected null. By default it is set equal to 1.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_selection_storey(Sx, Sy)
#' d_selection_storey(Sx, Sy, S=3)
#' d_selection_storey(Sx, Sy, S=c(3, 7:13))
#'
d_selection_storey = function(S_X, S_Y, S=NULL, alpha=0.1, lambda = 0.5, pvalue_only=FALSE){

  n = length(S_Y)
  m = length(S_X)
  pval_unsorted = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  r = order(pval_unsorted)
  pval = pval_unsorted[r]

  if(!pvalue_only){
   # Find h

   # For each level in the closed testing procedure, the worst case scenario Simes p-value is computed,
   # i.e., at level i in the closed testing procedure the Simes p-value is computed using
   # the n-i+1 greater conformal p-values.

   pval_simes = sapply(1:n, function(i)
     c(
       min(pval[i:n]/seq.int(from=1, to=n-i+1, by=1)),
       (1+sum(pval[i:n]>lambda))/((n-i+1)*(1-lambda))
     )
   )

   # In the last two steps of closed testing procedure
   # (the ones involving only elementary hypotheses or the hypotheses that are intersection of two elementary hypotheses)
   # Storey's estimator for the number of true null hypotheses is set equal to 1,
   # since it might take on values greater than 1.
   pval_simes[2,(n-1):n] <- 1

   d = sum(cumsum(pval_simes[1,] <= alpha/((n:1)*pval_simes[2,])) == 1:n)
   h = n - d

   # find d_S
   if (!is.null(S)){

     pval_S <- sort(pval_unsorted[S])
     s = length(pval_S)
     pval_notS <- sort(pval_unsorted[-S], decreasing = T)
     pvals = c(pval_notS, pval_S)
     o <- order(pvals)

     p_storey <- function(x) { ifelse( length(x) < 3, 1, (sum(x>lambda)+1)/(length(x)*(1-lambda)) ) * min(x*length(x)/(1:length(x)))  }

     d_S = 0
     for (i in 1:s){
       for (j in max(0,(h-s+i-1)):0 ){
         select = rep(F,length=n)
         select[c(0:j,((n-s)+i):n)] <- T
         notrejected <- p_storey( pval[select[o]] ) > alpha
         if (notrejected) { break }
       }
       if (notrejected) { break }
       d_S = d_S + 1
     }

   } else { d_S = d }


   # Compute global p-value and p-value for the selected null
   hom = hommel::hommel(pval_unsorted)
   # Compute p-value for the global null
   pval.global = hommel::localtest(hom)

   if(is.null(S)){
     # Compute p-value for the selected null
     pval.selection = pval.global

   } else {
     # Compute p-value for the selected null
     pval.selection = hommel::localtest(hom, ix=S)
   }

  } else {
    # Compute global p-value and p-value for the selected null
    hom = hommel::hommel(pval_unsorted)
    # Compute p-value for the global null
    pval.global = hommel::localtest(hom)
 }

  out = list("lower.bound" = d_S, 
             "global.p.value" = pval.global, 
             "S"=S, 
             "selection.p.value" = pval.selection)

  return(out)

}
