## Function to make density function monotone increasing
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
  g.mon <- function(x) g.hat.u(x) / integral
  
  return(g.mon)
}

## Function to make density function monotone decreasing
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
  g.mon <- function(x) g.hat.u(x) / integral
  
  return(g.mon)
}

## Function to choose between increasing and decreasing monotonicity
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

## Function to create CDF using interpolation
create_cdf_interpolated <- function(pdf_func, lower_bound=0, upper_bound=1, grid_points = 1000) {
    ## Generate a sequence of x values within the specified bounds
    x_grid <- seq(lower_bound, upper_bound, length.out = grid_points)
    
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

## Inverse transform sampling function
inverse_transform_sampling <- function(cdf_function, n_samples=1000, lower_bound=0, upper_bound=1) {
  ## Inverse transform sampling function
  iCDF <- GoFKernel::inverse(cdf_function, 0, 1)

  ## Generate random samples through inverse CDF transformation
  u_samples <- runif(n_samples)
  samples <- sapply(u_samples, function(u) iCDF(u))
  return(samples)
}

kernel_smoothed_pdf_function <- function(cdf_function, n_samples=1000) {
  ## Generate random samples using the inverse transform sampling
  samples <- inverse_transform_sampling(cdf_function)

  ## Kernel Density Estimation on the sampled values
  density_est <- density(samples, bw=0.05, kernel="gaussian", from=0, to=1)

  ## Normalize the PDF
  integral_value <- sum(density_est$y) * diff(density_est$x[1:2])
  normalized_pdf_values <- density_est$y / integral_value

  ## Create a function to return the normalized PDF
  normalized_pdf_function <- function(x) {
    approx(density_est$x, normalized_pdf_values, xout = x)$y
  }

  return(list(pdf_function = normalized_pdf_function))
}


fit_beta_mixture <- function(data, num_starts = 10, monotone=FALSE) {

    ## Log-likelihood function
    log_likelihood <- function(params, data) {
        ## Extract parameters
        alpha <- params[1]
        beta <- params[2]
        lambda <- params[3]

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

estimate_g <- function(scores_reference, scores_pooled, method="betamix", monotone=FALSE) {

    ## Transform the reference scores to make them approximately uniform
    null.fit <- fitdistrplus::fitdist(scores_reference, "beta", start = list(shape1 = 0.999, shape2 = 0.999))
    F.hat <- function(x) pbeta(x, null.fit$estimate[[1]], null.fit$estimate[[2]])
    scores_pooled = F.hat(scores_pooled)

    monotonicity <- NULL
    
    ## Fit beta mixture model
    betafit <- fit_beta_mixture(scores_pooled)
    g.hat <- function(x) dbeta(x, betafit$shape1, betafit$shape2)
    CDF.hat <- function(x) pbeta(x, betafit$shape1, betafit$shape2)
    if(method=="betamix") {
        if(monotone) {
            ## Decide between monotone increasing and decreasing approximation
            tmp <- choose_best_monotonic_density(g.hat)
            g.hat <- tmp$density
            monotonicity <- tmp$monotonicity
            ## Compute CDF
            CDF.hat <- create_cdf_interpolated(g.hat)
        }
    ## } else if (method=="mixmodel") {

    ##     ## Fit mixture model and evaluate CDF values on a grid
    ##     model <- mixmodel::mix.model(scores_pooled, method = "fixed" , c.n = .05*log(log(length(scores_pooled))), gridsize = 600)
    ##     ##model = hush(suppressWarnings(mixmodel::cv.mix.model(scores_pooled, cn.length=50)))

    ##     if(is.null(monotone)) {
    ##         CDF.values <- model$Fs.hat

    ##         ## Define CDF function via monotone increasing interpolation
    ##         CDF.values$y[1] = 0
    ##         CDF.values$y[length(CDF.values$y)] = 1
    ##         CDF.hat <- splinefun(c(0, CDF.values$x, 1), c(0, CDF.values$y, 1), method = "monoH.FC")

    ##         ## # DEBUG
    ##         ## x.grid <- seq(0,1,length.out=1000)
    ##         ## plot(x.grid, CDF.hat(x.grid), "l")

    ##         ## Estimate PDF from CDF via kernel smoothing
    ##         g.hat <- kernel_smoothed_pdf_function(CDF.hat)$pdf_function
    ##     } else if (monotone=="increasing") {
    ##         g.hat.values <- mixmodel::den.mix.model(model, dec.density=FALSE)
    ##         plot(g.hat.values$x, g.hat.values$y)
    ##         g.hat.u <- splinefun(g.hat.values$x[-c(1,length(g.hat.values$x))], g.hat.values$y[-c(1,length(g.hat.values$y))], method = "monoH.FC")
    ##         ## Normalize the density
    ##         integral <- integrate(g.hat.u, 0, 1)$value
    ##         g.hat <- function(x) g.hat.u(x) / integral
    ##         CDF.hat <- NULL
    ##     } else if (monotone=="decreasing") {
    ##         g.hat.values <- mixmodel::den.mix.model(model, dec.density=TRUE)
    ##         g.hat.u <- splinefun(g.hat.values$x[-c(1,length(g.hat.values$x))], g.hat.values$y[-c(1,length(g.hat.values$y))], method = "monoH.FC")
    ##         ## Normalize the density
    ##         integral <- integrate(g.hat.u, 0, 1)$value
    ##         g.hat <- function(x) g.hat.u(x) / integral
    ##         CDF.hat <- NULL
    ##     } else {
    ##         print("Error! Unknown parameter.")
    ##         g.hat <- NULL
    ##     }
    } else {
        stop("Error: unknown method!")
        g.hat <- NULL
        CDF.hat <- NULL
    }


  return(list(pdf = g.hat, cdf = CDF.hat, monotonicity=monotonicity))
}

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

compute.global.pvalue.shirashi.adaptive <- function(S_X, S_Y, prop_cal=0.5, num_mc=1000, fit.method="betamix", monotone=FALSE) {

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
    g.hat <- estimate_g(S_X1, S_pooled, method=fit.method, monotone=monotone)$pdf

    ## Carry out test by comparing S_X2 to S_Y
    p.val <- compute.global.pvalue.shirashi(S_X2, S_Y, g.hat, num_mc=num_mc)

    return(p.val)
}
