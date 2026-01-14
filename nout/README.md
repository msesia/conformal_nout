
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nout

<!-- badges: start -->
<!-- badges: end -->

The R package ***nout*** implements the closed testing shortcuts used by
ACODE.

## Installation

Install the package locally from this folder:

``` r
install.packages("devtools")  # if needed
devtools::install(".")
```

## Example

This is a basic example on synthetic data, where inlier scores are
generated from a Uniform(0,1) and outlier scores are generated under
Lehmann’s alternative of order k=2. All observations in the test set are
outliers.

``` r
library(nout)

set.seed(321)

# Auxiliary functions to generate outliers scores
g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
rg2 = function(rnull, k=2) max(rnull(k))

# Generating calibration and test score vectors
# Outlier distribution is the Lehmann's alternative of order k=2
X = runif(20)
Y = replicate(20, rg2(rnull=runif))
```

### Lower bound without selection

``` r
# Confidence lower-bound for the number of outliers using Closed Testing with local Simes test
res1 = d_selection_simes(X, Y,alpha = 0.1)
res1
#> $lower.bound
#> [1] 0
#> 
#> $global.p.value
#> [1] 0.3174603
#> 
#> $S
#> NULL
#> 
#> $selection.p.value
#> [1] 0.3174603

# Confidence lower-bound for the number of outliers using Closed Testing with local Fisher's method
res2 = d_selection_fisher(X, Y, B=100, n_perm=0, alpha=0.1)
res2
#> $lower.bound
#> [1] 2
#> 
#> $global.pvalue
#> [1] 0.05081268
#> 
#> $S
#> NULL
#> 
#> $selection.p.value
#> [1] 1

# Confidence lower-bound for the number of outliers using Closed Testing with local Wilcoxon sum-rank test 
res3 = d_selection_higher(X, Y, local_test="WMW", n_perm=0, B=100, alpha=0.1)
res3
#> $lower.bound
#> [1] 0
#> 
#> $global.pvalue
#> [1] 0.00549985
#> 
#> $S
#> NULL
#> 
#> $selection.p.value
#> [1] 1

# Confidence lower-bound for the number of outliers using Closed Testing with local higher order
# Wilcoxon sum-rank test of order k=2
res4 = d_selection_higher(X, Y, local_test="higher", k=2, n_perm=0, B=100, alpha=0.1)
res4
#> $lower.bound
#> [1] 7
#> 
#> $global.pvalue
#> [1] 0.005205133
#> 
#> $S
#> NULL
#> 
#> $selection.p.value
#> [1] 1

# Confidence lower-bound for the number of outliers using Closed Testing with local oracle Shiraishi test
res5 = d_selection_G(X, Y, g.hat = g2, monotone=TRUE, B=100, alpha=0.1)
res5
#> $lower.bound
#> [1] 7
#> 
#> $global.pvalue
#> [1] 0.005510923
#> 
#> $S
#> NULL
#> 
#> $selection.p.value
#> [1] 1
```

### Lower bound with selection

The selection set is randomly chosen among the test index set.

``` r
S = sample(1:20, size = 16)

# Confidence lower-bound for the number of outliers using Closed Testing with local Simes test
resS1 = d_selection_simes(X, Y, S=S, alpha = 0.1)
resS1
#> $lower.bound
#> [1] 0
#> 
#> $global.p.value
#> [1] 0.3174603
#> 
#> $S
#>  [1]  6  3  8 11  1 20 15 10 13 14  2  5 16 12  7 17
#> 
#> $selection.p.value
#> [1] 0.3373016

# Confidence lower-bound for the number of outliers using Closed Testing with local Fisher's method
resS2 = d_selection_fisher(X, Y, S=S, B=100, n_perm=0, alpha=0.1)
resS2
#> $lower.bound
#> [1] 1
#> 
#> $global.pvalue
#> [1] 0.05081268
#> 
#> $S
#>  [1]  6  3  8 11  1 20 15 10 13 14  2  5 16 12  7 17
#> 
#> $selection.p.value
#> [1] 1

# Confidence lower-bound for the number of outliers using Closed Testing with local Wilcoxon sum-rank test 
resS3 = d_selection_higher(X, Y, S=S, local_test="WMW", n_perm=0, B=100, alpha=0.1)
resS3
#> $lower.bound
#> [1] 5
#> 
#> $global.pvalue
#> [1] 0.006510641
#> 
#> $S
#>  [1]  6  3  8 11  1 20 15 10 13 14  2  5 16 12  7 17
#> 
#> $selection.p.value
#> [1] 1

# Confidence lower-bound for the number of outliers using Closed Testing with local higher order
# Wilcoxon sum-rank test of order k=2
resS4 = d_selection_higher(X, Y, S=S, local_test="higher", k=2, n_perm=0, B=100, alpha=0.1)
resS4
#> $lower.bound
#> [1] 5
#> 
#> $global.pvalue
#> [1] 0.04111738
#> 
#> $S
#>  [1]  6  3  8 11  1 20 15 10 13 14  2  5 16 12  7 17
#> 
#> $selection.p.value
#> [1] 1

# Confidence lower-bound for the number of outliers using Closed Testing with local oracle Shiraishi test
resS5 = d_selection_G(X, Y, S=S, g.hat = g2, monotone=TRUE, B=100, alpha=0.1)
resS5
#> $lower.bound
#> [1] 5
#> 
#> $global.pvalue
#> [1] 0.005483834
#> 
#> $S
#>  [1]  6  3  8 11  1 20 15 10 13 14  2  5 16 12  7 17
#> 
#> $selection.p.value
#> [1] 1
```
