rm(list=ls())

results_file <- "results/tabA2/simes_perm_results.csv"
dir.create("results/tabA2", showWarnings = FALSE, recursive = TRUE)

require(combinat)

MS <- c(9,14,19,24,29,34,39,44,49,54)
SIZE_SIMES_PERM <- vector()
SIZE_SIMES <- vector()
CRIT_VAL <- vector()

alpha = 0.1 
n = 3

for (i in 1:length(MS)){
    
    m = MS[i]
    N = m + n
    Z <- 1:N
    perms = combn(1:N,m)
    B = ncol(perms)
    stat<-vector()
    for (b in 1:B){
        Zperm = c(Z[perms[,b]], Z[-perms[,b]])
        p_conformal = sapply(1:n, function(i)
            wilcox.test(x=Zperm[m+i],y=Zperm[1:m],
                        exact=T,
                        alternative="greater")$p.value)
        stat[b] = min( sort(p_conformal)*n/(1:n) )
    }
    ord_stat <- sort(unique(stat))
    q1 <- as.numeric(quantile(stat, probs=alpha, type=1))
    
    if ( mean(stat <= q1) > alpha ) {
        
        k <- which.min(ord_stat < q1)
        crit <- ord_stat[k-1]
        
    } else { crit = q1 }
    
    size <- mean(stat <= crit)
    SIZE_SIMES_PERM [i] <- size
    CRIT_VAL[i] <- crit
    SIZE_SIMES[i] <- ecdf(stat)(alpha)
}

out_df <- data.frame(
    m               = MS,
    size_simes      = SIZE_SIMES,
    size_simes_perm = SIZE_SIMES_PERM,
    crit_val        = CRIT_VAL
)
write.csv(out_df, file = results_file, row.names = FALSE)
