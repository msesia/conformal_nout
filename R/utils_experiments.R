## Function to compute global p-values using various statistical tests
##
## Args:
##   data: List. A list containing:
##     - scores.cal: Numeric vector. Calibration set scores.
##     - scores.test: Numeric vector. Test set scores.
##     - outlier.test: Numeric vector. Indicator of outliers in the test set (0 = inlier, 1 = outlier).
##   alternative: String. The type of distribution for outliers in the test set.
##
## Returns:
##   A data frame with the computed global p-values for each method.
run_global_testing <- function(data, alternative=NULL) {

    S_X = data$scores.cal
    S_Y = data$scores.test

    ## Compute global p-value with Fisher's test
    pval.fisher <- d_selection_fisher(S_X = S_X, S_Y = S_Y, n_perm = 0, pvalue_only = TRUE)$global.pvalue

    ## Compute global p-value with WMW test
    pval.wmw <- d_selection_higher(S_X = S_X, S_Y = S_Y, local.test = "WMW", n_perm = 0, pvalue_only = TRUE)$global.pvalue

    ## Compute global p-value with WMW test (k=2)
    pval.wmw.2 <- d_selection_higher(S_X = S_X, S_Y = S_Y, local.test = "higher", k = 2, n_perm = 0, pvalue_only = TRUE)$global.pvalue

    ## Compute global p-value with WMW test (k=3)
    pval.wmw.3 <- d_selection_higher(S_X = S_X, S_Y = S_Y, local.test = "higher", k = 3, n_perm = 0, pvalue_only = TRUE)$global.pvalue

    ## Compute global p-value with Shirashi's approach (using oracle density)
    if(!is.null(alternative)) {
        density_oracle <- function(x) density_scores(x, alternative)
        pval.g.oracle <- compute.global.pvalue.shirashi(S_X, S_Y, density_oracle, num_mc=1000)
    } else {
        pval.g.oracle <- NA
    }

    ## Apply Shirashi's method using g-hat estimated through beta mixture
    pval.g.hat.1 <- compute.global.pvalue.shirashi.adaptive(S_X, S_Y, prop_cal=0.5, num_mc=1000, fit.method="betamix")

        ## Apply Shirashi's method using g-hat estimated through beta mixture (increasing)
    pval.g.hat.2 <- compute.global.pvalue.shirashi.adaptive(S_X, S_Y, prop_cal=0.5, num_mc=1000, fit.method="betamix", monotone="increasing")

    ## Apply Shirashi's method using g-hat estimated through mixmodel
    pval.g.hat.3 <- compute.global.pvalue.shirashi.adaptive(S_X, S_Y, prop_cal=0.5, num_mc=1000, fit.method="mixmodel")

    ## Apply Shirashi's method using g-hat estimated through mixmodel (increasing)
    pval.g.hat.4 <- compute.global.pvalue.shirashi.adaptive(S_X, S_Y, prop_cal=0.5, num_mc=1000, fit.method="mixmodel", monotone="increasing")

    ## Create a data frame with the p-values
    pval_df <- tibble::tibble(
                           Method = c("Fisher", "WMW", "WMW_k2", "WMW_k3", "Shirashi_oracle", "Shirashi_ghat_betamix", "Shirashi_ghat_betamix_inc",
                                      "Shirashi_ghat_mixmodel", "Shirashi_ghat_mixmodel_inc"),
                           p.value = c(pval.fisher, pval.wmw, pval.wmw.2, pval.wmw.3, pval.g.oracle, pval.g.hat.1, pval.g.hat.2, pval.g.hat.3, pval.g.hat.4)
                       )

    return(pval_df)
}


run_outlier_enumeration <- function(data, alternative=NULL) {

    S_X = data$scores.cal
    S_Y = data$scores.test

    ## Estimate the number of outliers with Fisher's test
    res.fisher <- d_selection_fisher(S_X = S_X, S_Y = S_Y, n_perm = 0)

    ## Compute global p-value with WMW test
    res.wmw <- d_selection_higher(S_X = S_X, S_Y = S_Y, local.test = "WMW", n_perm = 0)

    ## Compute global p-value with WMW test (higher order)
    res.wmw.2 <- d_selection_higher(S_X = S_X, S_Y = S_Y, local.test = "higher", k=2, n_perm = 0)

    ## Compute global p-value with WMW test (higher order)
    res.wmw.3 <- d_selection_higher(S_X = S_X, S_Y = S_Y, local.test = "higher", k=3, n_perm = 0)

    ## Compute global p-value with Shirashi's approach (using oracle density)
    if(!is.null(alternative)) {
        density_oracle <- function(x) density_scores(x, alternative)
        ## Make sure the density is monotone increasing
        density_oracle_m <- make_density_monotone(density_oracle)       
        res.g.oracle <- d_selection_G2(S_X, S_Y, g.oracle=density_oracle_m, monotonicity="increasing", alpha=0.1, n_perm=0, B=10^3, B_MC=10^3, seed=123)
        d.g.oracle <- res.g.oracle$lower.bound
        pval.g.oracle <- res.g.oracle$p.value
    } else {
        d.g.oracle <- NA
        pval.g.oracle <- NA
    }

    ## Apply Shirashi's method using g-hat estimated through beta mixture (increasing)
    res.g.hat.1 <- d_selection_G2(S_X, S_Y, g.oracle=NULL, fit.method="betamix", monotonicity="increasing", prop.cal=0.5, alpha=0.1, n_perm=0, B=10^3, B_MC=10^3, seed=123)

    ## Apply Shirashi's method using g-hat estimated through mixmodel (increasing)
    res.g.hat.2 <- d_selection_G2(S_X, S_Y, g.oracle=NULL, fit.method="mixmodel", monotonicity="increasing", prop.cal=0.5, alpha=0.1, n_perm=0, B=10^3, B_MC=10^3, seed=123)

    ## Create a data frame with the p-values
    df <- tibble::tibble(
                           Method = c("Fisher", "WMW", "WMW_k2", "WMW_k3", "Shirashi_oracle_inc", "Shirashi_g_hat_betamix_inc", "Shirashi_g_hat_mixmodel_inc"),
                           Lower = c(res.fisher$lower.bound, res.wmw$lower.bound, res.wmw.2$lower.bound, res.wmw.3$lower.bound,
                                     d.g.oracle, res.g.hat.1$lower.bound, res.g.hat.2$lower.bound),
                           p.value = c(res.fisher$p.value, res.wmw$p.value, res.wmw.2$p.value, res.wmw.3$p.value,
                                       pval.g.oracle, res.g.hat.1$p.value, res.g.hat.2$p.value)
                       )
    return(df)
}
