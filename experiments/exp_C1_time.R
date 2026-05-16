library(ggplot2)

results_file <- "results/figA1/time_results.csv"
dir.create("results/figA1", showWarnings = FALSE, recursive = TRUE)

require(hommel)
require(nout)

set.seed(123)

size = c(5000,10000,15000,20000,25000)
TOTAL_TIME = matrix(NA,nrow=length(size),ncol=6)

for (kk in 1:length(size) ){
    
    m = n = size[kk]
    N = n+m
    theta = 0.5
    n0 = n*(1-theta)
    alpha = 0.1
    lambda = 0.5
    mu = sqrt(2*log(n))
    g = function(x, mu=sqrt(2*log(n)) ) ifelse(x<1 & x>0, dnorm( qnorm(x) , mean = mu )*(1/dnorm(qnorm(x)) ), 0)
    crit_fisher = sapply(1:n, function(h) sqrt(1+(h/m)) * stats::qchisq(p=1-alpha, df=2*h) - 2 * (sqrt(1+(h/m))-1) * h )
    crit_wmw = sapply(1:n, function(h) as.double(stats::qnorm(alpha, mean=m*h/2, sd = sqrt(m*h*(m+h+1)/12), lower.tail = F)))
    
    B = 3
    TIME <- matrix(NA,nrow=B,ncol=6)
    
    for (b in 1:B){
        
        t <- system.time({
            X = sort(rnorm(m))
            Y = sort(c(rnorm(n0), rnorm(n-n0, mean=mu)))
            W = c(X,Y)
            R = rank(W, ties.method="first")
            rY <- R[(m+1):(m+n)]
            t  <- 1:n - 1L
            r_tilde <- rY - t
            pval <- rev( (m + 2 - r_tilde) / (m + 1) )
        })
        
                                        # simes
        t_simes <- system.time({
            hom = hommel::hommel(pval)
            d = hommel::discoveries(hom, alpha = alpha)
        })
                                        # storey
        t_storey <- system.time({
            
            gt <- integer(n + 2L)
            for (k in n:1L) gt[k] <- gt[k + 1L] + (pval[k] > lambda)
            
            h <- n
            while (h >= 1L) {
                start <- n - h + 1L
                min_ratio <- Inf
                idx <- start
                for (j in 1L:h) {
                    r <- pval[idx] / j
                    if (r < min_ratio) min_ratio <- r
                    idx <- idx + 1L
                }
                lhs <- h * min_ratio
                rhs <- alpha * (h * (1 - lambda)) / (1 + gt[start])
                if (lhs <= rhs) {
                    h <- h - 1L
                } else {
                    break
                }
            }
            
        })
        
                                        # shiraishi
        t_shiraishi <- system.time({
            d = 0
            for (i in n:1){
                a_i = g( (1:(m+i)) / ((m+i)+1) )
                stat_i = sum(  a_i[ R[(m+1):(m+i)] ]  )
                mu_i = i*mean(a_i)
                var_i = i*m*sum((a_i-mean(a_i))^2)/((m+i)*((m+i)-1))
                d = d  + ( stat_i > mu_i + qnorm(1-alpha)*sqrt(var_i) )
                if (d < n - i + 1){ break }
            }
        })
                                        # fisher
        t_fisher <- system.time({
            sumSome::discoveries(sumSome::sumStatsPar(g = -2 * log(pval), alpha = alpha, 
                                                      cvs = crit_fisher))
        })
                                        # wmw
        t_wmw <- system.time({
            sumSome::discoveries(sumSome::sumStatsPar(g = rY - seq_len(n), alpha = alpha, cvs = crit_wmw))
        })
        
        TIME[b,1] = t["elapsed"]
        TIME[b,2] = t["elapsed"] + t_simes["elapsed"]
        TIME[b,3] = t["elapsed"] + t_storey["elapsed"]
        TIME[b,4] = t["elapsed"] + t_shiraishi["elapsed"]
        TIME[b,5] = t["elapsed"] + t_fisher["elapsed"]
        TIME[b,6] = t["elapsed"] + t_wmw["elapsed"]
    }
    print(kk)
    TOTAL_TIME[kk,] <- apply(TIME,2,median)
    
}

TOTAL_TIME_MS <- TOTAL_TIME*1000

out_df <- data.frame(
    n = size,
    Preparatory_step = TOTAL_TIME_MS[,1],
    Simes            = TOTAL_TIME_MS[,2],
    Storey           = TOTAL_TIME_MS[,3],
    Shiraishi        = TOTAL_TIME_MS[,4],
    Fisher           = TOTAL_TIME_MS[,5],
    WMW              = TOTAL_TIME_MS[,6]
)
write.csv(out_df, file = results_file, row.names = FALSE)
