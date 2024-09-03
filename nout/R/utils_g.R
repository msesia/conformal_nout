# --------------------------------------------------------------- #
#  Implementation of Shiraishi test with asymptotic distribution  #
# --------------------------------------------------------------- #


#' standardize_to_uniform
#' 
#' @description
#' It rank the data and scale ranks to [0, 1].
#' 
#' @param X1 : a vector
#' @param X2 : a vector
#'
#' @return A list of two vectors corresponding to the scaled ranks of the two input vectors
#'
standardize_to_uniform <- function(X1, X2) {
  # Rank the data and scale ranks to [0, 1]
  n1 <- length(X1)
  X.pool <- c(X1, X2)
  n.pool <- length(X.pool)
  uniform_data <- rank(X.pool, ties.method = "average") / (n.pool + 1)
  X1.s <- uniform_data[1:n1]
  X2.s <- uniform_data[(n1+1):n.pool]
  return(list(X1=X1.s, X2=X2.s))
}


#' make_density_monotone_increasing
#'
#' @description
#' It makes a density function monotone increasing using Isotonic Regression
#' 
#' @param g : a density function
#'
#' @return A monotone increasing density function that is the output of Isotonic Regression
#'
make_density_monotone_increasing <- function(g) {
  tol = 1e-3
  x.grid <- seq(tol,1-tol,length.out=1000)

  ## Obtain density estimates from the model
  g.hat.values <- list(x=x.grid, y=g(x.grid))

  ## Create a monotone interpolation function
  iso_fit <- isoreg(g.hat.values$x, g.hat.values$y)
  g.hat.u <- approxfun(iso_fit$x, iso_fit$yf, method = "linear", rule = 2)

  ## Normalize the density
  integral <- integrate(g.hat.u, 0, 1, stop.on.error=FALSE)$value
  g.mon <- function(x) g.hat.u(x) / pmax(1e-6, integral)

  return(g.mon)
}

#' make_density_monotone_decreasing
#'
#' @description
#' It makes a density function monotone decreasing using Isotonic Regression
#' 
#' @param g : a density function
#'
#' @return A monotone decreasing density function that is the output of Isotonic Regression
#'
make_density_monotone_decreasing <- function(g) {
  tol = 1e-3
  x.grid <- seq(tol, 1 - tol, length.out = 1000)

  ## Obtain density estimates from the model
  g.hat.values <- list(x = x.grid, y = g(x.grid))

  ## Create a decreasing monotone interpolation function
  # Reverse the order of y values and apply isoreg to ensure monotonicity
  reverse_y <- rev(g.hat.values$y)
  iso_fit <- isoreg(g.hat.values$x, reverse_y)
  decreasing_y <- rev(iso_fit$yf)

  # Create an interpolation function from the decreasing values
  g.hat.u <- approxfun(g.hat.values$x, decreasing_y, method = "linear", rule = 2)

  ## Normalize the density
  integral <- integrate(g.hat.u, 0, 1, stop.on.error=FALSE)$value
  g.mon <- function(x) g.hat.u(x) / pmax(1e-6, integral)

  return(g.mon)
}



#' choose_best_monotonic_density
#'
#' @description
#' It chooses between increasing and decreasing monotonicity 
#' 
#' @param g : a density function
#'
#' @return A density function that is the output of Isotonic Regression and 
#' is monotone decreasing or increasing based on the lower Residual Sum of Squares
#'
choose_best_monotonic_density <- function(g) {
  tol = 1e-3
  x.grid <- seq(tol, 1 - tol, length.out = 1000)

  # Generate monotone increasing and decreasing density functions
  g.hat.inc <- make_density_monotone_increasing(g)
  g.hat.dec <- make_density_monotone_decreasing(g)

  # Calculate residual sum of squares (RSS) for increasing and decreasing densities
  rss.inc <- sum((g.hat.inc(x.grid) - g(x.grid))^2)
  rss.dec <- sum((g.hat.dec(x.grid) - g(x.grid))^2)

  # Decide based on lower RSS
  if (rss.inc < rss.dec) {
    list(density = g.hat.inc, monotonicity = "increasing")
  } else {
    list(density = g.hat.dec, monotonicity = "decreasing")
  }
}


#' create_cdf_interpolated
#'
#' @description
#' It creates cumulative density function using linear interpolation within a grid of points 
#' 
#' @param pdf_func : density function
#' @param lower_bound : lower bound of the density support 
#' @param upper_bound : upper bound of the density support 
#' @param grid_points : number of points where to compute the cumulative density function via integration of the density function
#'
#' @return A function which is the interpolated cumulative density function
#'
create_cdf_interpolated <- function(pdf_func, lower_bound=0, upper_bound=1, grid_points = 1000) {
    ## Generate a sequence of x values within the specified bounds
    tol = 1e-4
    x_grid <- seq(lower_bound+tol, upper_bound-tol, length.out = grid_points)

    ## Compute CDF values on the grid using numerical integration
    cdf_values <- numeric(grid_points)
    for (i in 1:grid_points) {
        cdf_values[i] <- integrate(pdf_func, lower_bound, x_grid[i], stop.on.error=FALSE)$value
    }

    ## Create interpolation function
    cdf_func <- suppressWarnings(approxfun(c(lower_bound,x_grid,upper_bound), c(0,cdf_values,1), rule = 2))

    return(cdf_func)
}


hush <- function(code){
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
}


#' inverse_transform_sampling
#'
#' @description
#' Given a cumulative density function, it generates a specific number of samples 
#' from the specified distribution via inverse transform sampling
#' 
#' @param pdf_func : cumulative density function
#' @param n_samples : number of samples to be generated
#' @param lower_bound : lower bound of the density support 
#' @param upper_bound : upper bound of the density support 
#'
#' @return A vector of samples from the desired distribution of the prespecified length
#'
inverse_transform_sampling <- function(cdf_function, n_samples=1000, lower_bound=0, upper_bound=1) {
  ## Inverse transform sampling function
  iCDF <- GoFKernel::inverse(cdf_function, 0, 1)

  ## Generate random samples through inverse CDF transformation
  u_samples <- runif(n_samples)
  samples <- sapply(u_samples, function(u) iCDF(u))
  return(samples)
}




#' kernel_smoothed_pdf_function
#' 
#' @description
#' Given the cumulative density function, it returns a normalized Kernel Density Estimate of the density function
#' 
#' @param cdf_function : a cumulative density function
#' @param n_samples : number of samples to be generated to estimate the density function via Kernel Density Estimation
#'
#' @return A function which is the normalized Kernel Density Estimate of the density function
#' 
kernel_smoothed_pdf_function <- function(cdf_function, n_samples=1000) {
  ## Generate random samples using the inverse transform sampling
  samples <- inverse_transform_sampling(cdf_function)

  ## Kernel Density Estimation on the sampled values
  density_est <- density(samples, bw=0.05, kernel="gaussian", from=0, to=1)

  ## Normalize the PDF
  integral_value <- sum(density_est$y) * diff(density_est$x[1:2])
  normalized_pdf_values <- density_est$y / pmax(1e-6, integral_value)

  ## Create a function to return the normalized PDF
  normalized_pdf_function <- function(x) {
    approx(density_est$x, normalized_pdf_values, xout = x)$y
  }

  return(list(pdf_function = normalized_pdf_function))
}




#' fit_mixmodel
#'
#' @description
#' It estimates the CDF from a vector of unidimensional scores and then, from the estimated CDF, it returns the PDF estimated via kernel smoothing
#' 
#' @param data : vector of univariate scores
#'
#' @return A function which the PDF estimated via kernel smoothing from the score vector in input
#' 
fit_mixmodel <- function(data) {

    ## Fit mixture model and evaluate CDF values on a grid
    model <- mixmodel::mix.model(data, method = "fixed" , c.n = .05*log(log(length(data))), gridsize = 600)
    ##model = hush(suppressWarnings(mixmodel::cv.mix.model(data, cn.length=50)))

    CDF.values <- model$Fs.hat

    CDF.values$y[1] = 0
    CDF.values$y[length(CDF.values$y)] = 1
    CDF.hat <- splinefun(c(0, CDF.values$x, 1), c(0, CDF.values$y, 1), method = "monoH.FC")

    ## Estimate PDF from CDF via kernel smoothing
    g.hat <- kernel_smoothed_pdf_function(CDF.hat)$pdf_function
    return(g.hat)
}




#' fit_beta_mixture
#' 
#' @description
#' It estimates the parameters of a mixture model of two distributions, 
#' where the first one is a Standard Uniform and the second one is a Beta distribution 
#' whose parameters are both unknown and also the mixing proportion is unknown. 
#' All the parameters are optimized minimizing the negative log-likelihood,
#' allowing for box constraints.
#' 
#' @param data : score vector
#' @param num_starts : number of iterations for parameter optimization
#'
#' @return A list of three numbers that are, respectively, the estimates of 
#' the first parameter (alpha) of a Beta distribution, of the second parameter (beta) of a Beta distribution and
#' of the mixing proportion (lambda) of the mixture distribution.
#' 
fit_beta_mixture <- function(data, num_starts = 10) {

    ## Log-likelihood function
    log_likelihood <- function(params, data) {
        ## Extract parameters
        alpha <- params[1] # first parameter of Beta distribution
        beta <- params[2] # first parameter of Beta distribution
        lambda <- params[3] # mixing proportion

        ## Ensure parameters are in valid range
        if (alpha <= 0 || beta <= 0 || lambda <= 0 || lambda >= 1) return(-Inf)

        ## Density of Beta distribution
        beta_density <- dbeta(data, shape1 = alpha, shape2 = beta)

        ## Density of Uniform distribution
        uniform_density <- dunif(data, min = 0, max = 1)

        ## Mixture density
        mixture_density <- lambda * uniform_density + (1 - lambda) * beta_density

        ## Log-likelihood
        log_likelihood <- sum(log(mixture_density))

        return(-log_likelihood) # Return negative log-likelihood for minimization
    }

    ## Make sure there are no 0's or 1's
    tol <- 1e-3
    data <- pmax(data, tol)
    data <- pmin(data, 1 - tol)

    ## Set up multiple starts for initial parameters
    best_fit <- NULL
    best_log_likelihood <- Inf

    lower <- c(0.01, 0.01, 0.01)
    upper <- c(Inf, Inf, 0.99)

    for (i in 1:num_starts) {
        ## Random initial parameter guesses
        initial_params <- c(runif(1, 0.5, 2), runif(1, 0.5, 2), runif(1, 0.01, 0.99))

        ## Optimization to find the best parameters
        fit <- optim(
            par = initial_params,
            fn = log_likelihood,
            data = data,
            method = "L-BFGS-B",
            lower = lower,
            upper = upper
        )

        ## Check if the current fit is the best one
        if (fit$value < best_log_likelihood) {
            best_log_likelihood <- fit$value
            best_fit <- fit
        }
    }

    ## Extract best fitted parameters
    alpha_hat <- best_fit$par[1]
    beta_hat <- best_fit$par[2]
    lambda_hat <- best_fit$par[3]

    return(list(shape1 = alpha_hat, shape2 = beta_hat, lambda = lambda_hat))
}




#' estimate_g
#'
#' @description
#' Given a two-component mixture model where the first distribution is a Standard Uniform 
#' and the second one is a Beta distribution, it estimates the parameter of the Beta distribution 
#' and returns the PDF and the CDF.
#' 
#' @param scores_reference : calibration score vector of inliers.
#' @param scores_pooled : pooled score vector of calibration and test scores.
#' @param method : character value indicating the method to be used to estimate the PDF. 
#' It can be either "betamix" or "mixmodel".
#' @param monotone : logical value. If \code{TRUE} the estimated PDF via the specified method 
#' is made monotone. The direction (increasing or decreasing) is automatically chosen from the data.
#' 
#' @return A list of three elements which are the estimated PDF and CDF and the monotonicity 
#' used in the estimation process.
#' 
estimate_g <- function(scores_reference, scores_pooled, method="betamix", monotone=FALSE) {
    ## Transform the reference scores to make them approximately uniform
    null.fit <- fitdistrplus::fitdist(as.numeric(scores_reference), "beta", start = list(shape1 = 0.999, shape2 = 0.999))
    F.hat <- function(x) pbeta(x, null.fit$estimate[[1]], null.fit$estimate[[2]])
    scores_pooled = F.hat(scores_pooled)

    monotonicity <- NULL

    if(method=="betamix") {
        ## Fit beta mixture model
        betafit <- fit_beta_mixture(scores_pooled)
        g.hat <- function(x) dbeta(x, betafit$shape1, betafit$shape2)
        g.hat.o <- g.hat
        CDF.hat <- function(x) pbeta(x, betafit$shape1, betafit$shape2)
    } else if (method=="mixmodel") {
        g.hat <- fit_mixmodel(scores_pooled)
    } else {
        stop("Error: unknown method!")
        g.hat <- NULL
        CDF.hat <- NULL
    }
    if(monotone) {
        ## Decide between monotone increasing and decreasing approximation
        tmp <- choose_best_monotonic_density(g.hat)
        g.hat <- tmp$density
        monotonicity <- tmp$monotonicity
    }
    ## Compute CDF
    CDF.hat <- create_cdf_interpolated(g.hat)

    return(list(pdf = g.hat, cdf = CDF.hat, monotonicity=monotonicity))
}




#' compute.global.pvalue.shirashi
#'
#' @description
#' Given the outlier density function and the calibration and test score vectors,
#' it computes the global *p*-value for the Shiraishi test.
#' 
#' @param S_X : calibration score vector.
#' @param S_Y : test score vector.
#' @param g : outlier density.
#' @param num_mc : number of Monte Carlo iterations.
#'
#' @return A number, the global *p*-value for the Shiraishi test.
#' 
compute.global.pvalue.shirashi <- function(S_X, S_Y, g, num_mc=1000) {
    m <- length(S_X)
    n <- length(S_Y)
    N <- m+n

    S_pool = c(S_X, S_Y)
    R = rank(S_pool)[(m+1):(m+n)]
    T.obs <- mean(sapply(1:num_mc, function(b) {
        U = sort(stats::runif(N))
        return(sum(g(U[R])))
    }))

    Z <- t(sapply(1:num_mc, function(b) {
        U = sort(stats::runif(N))
        return(g(U))
    }))
    Z.exp = colMeans(Z)

    mu <- mean(Z.exp)
    sigma <- sqrt(m*n*sum((Z.exp-mu)^2) / (N*(N-1)))

    if(sigma>0) {
        p.value = stats::pnorm(q=T.obs, mean=n*mu, sd = sigma, lower.tail = F)
    } else {
        p.value = ifelse(n*mu >= T.obs, runif(1), 0)
    }

    return(p.value)
}





#' compute.global.pvalue.shirashi
#'
#' @description
#' Given the calibration and test score vectors, it computes the global *p*-value 
#' for the Shiraishi test after estimating the outlier density function.
#' 
#' @param S_X : calibration score vector.
#' @param S_Y : test score vector.
#' @param prop_cal : proportion of inlier observations used for estimating the inlier distribution.
#' @param num_mc : number of Monte Carlo iterations.
#' @param fit_method : character value indicating the method for estimating the PDF. 
#' It can be either "betamix" or "mixmodel".
#' @param monotone : logical value. If \code{TRUE} the estimated PDF via the specified method 
#' is made monotone. The direction (increasing or decreasing) is automatically chosen from the data.
#'
#' @return A number, the global *p*-value for the Shiraishi test when the outlier density is estimated from the data.
#' 
compute.global.pvalue.shirashi.adaptive <- function(S_X, S_Y, prop_cal=0.5, num_mc=1000, fit_method="betamix", monotone=FALSE) {

    ## Split the reference scores and create pooled vector
    m = length(S_X)
    n = length(S_Y)
    m_2 = pmin(n, as.integer(round(prop_cal * m)))
    m_1 = m - m_2
    idx_X_1 = sample(m, m_1)
    idx_X_2 = setdiff(1:m, idx_X_1)
    S_X1 = S_X[idx_X_1]
    S_X2 = S_X[idx_X_2]
    S_pooled = sample(c(S_X2, S_Y))

    ## Estimate g-hat by comparing S_X1 to S_pooled
    g.hat <- estimate_g(S_X1, S_pooled, method=fit_method, monotone=monotone)$pdf

    ## Carry out test by comparing S_X2 to S_Y
    p.val <- compute.global.pvalue.shirashi(S_X2, S_Y, g.hat, num_mc=num_mc)

    return(p.val)
}
