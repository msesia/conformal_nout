## -----------------------------------------------------------------------------
## Required results : fig1
## Produces         : Figure 1, Figure A13, Table A4
## -----------------------------------------------------------------------------

source("utils_plotting_A.R")
library(ggh4x)

## ---- Plot function: 3-facet version (Median LB | Quantile.90 LB | Power) ----
make_plot_lower_bound_lhco <- function(fig.name, plot.n_train, plot.n_cal, tab.name="tab?", include.BH=FALSE,  plot.classifier="auto",
                                       plot.simple=FALSE, plot.naive=FALSE, save.plot=TRUE, alpha=0.1, n_test=2000, max_out=0.15) {
    if(plot.simple) {
        init_settings(idx.exclude=c(2,5,6))
    } else {
        init_settings(idx.exclude=c(6))
    }
    pow.str <- "Power (global null)"
    key.values <- c("Median", "Quantile.90", pow.str)
    key.labels <- c("Lower bound (median)", "Lower bound (90th quant.)", pow.str)
    disc.str <- "Discoveries (10% FDR)"
    df.lb <- results %>%
        filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier==plot.classifier) %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_higher_k2"),
                     names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method, tune_size, selection) %>%
        summarise(Median=median(lower_bound), Quantile.90=quantile(lower_bound, 0.9), SE=sd(lower_bound)/sqrt(n()), N=n()) %>%
        pivot_longer(c("Median", "Quantile.90"), names_to="Key", values_to="Value")
    if(plot.naive) {
        method.values.tmp <- c(method.values, "greedy")
        method.labels.tmp <- c(method.labels, "Cherry picking")
        color.scale.tmp <- c(color.by.label, "Cherry picking" = cbPalette[2])
        shape.scale.tmp <- c(shape.by.label, "Cherry picking" = 4)
        alpha.scale.tmp <- c(alpha.by.label, "Cherry picking" = 0.5)
        df.lb.greedy <- results %>%
            filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier!="auto") %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw", "lb_higher_k2"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, tune_size, selection, seed, Seed, Repetition) %>%
            summarize(lower_bound_greedy=max(lower_bound)) %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, tune_size, selection) %>%
            summarise(method="greedy", Median=median(lower_bound_greedy), Quantile.90=quantile(lower_bound_greedy, 0.9), SE=sd(lower_bound_greedy)/sqrt(n()), N=n()) %>%
            pivot_longer(c("Median", "Quantile.90"), names_to="Key", values_to="Value")
        df.lb <- df.lb %>%
            rbind(df.lb.greedy) %>%
            filter(method %in% method.values.tmp) %>%
            mutate(Method = factor(method, method.values.tmp, method.labels.tmp), Target="Closed testing")
    } else {
        method.values.tmp <- method.values
        method.labels.tmp <- method.labels
        color.scale.tmp <- color.by.label
        shape.scale.tmp <- shape.by.label
        alpha.scale.tmp <- alpha.by.label
        df.lb <- df.lb %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels), Target="Closed testing")
    }
    df.pow <- results %>%
        filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier==plot.classifier) %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_higher_k2"),
                     names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method, tune_size, selection) %>%
        summarise(Key=pow.str, Value=mean(lower_bound>0), SE=sd(lower_bound>0)/sqrt(n()), N=n())
    if(plot.naive) {
        df.pow.greedy <- results %>%
            filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier!="auto") %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw", "lb_higher_k2"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, tune_size, selection, seed, Seed, Repetition) %>%
            summarize(pow_greedy=max(lower_bound)) %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, tune_size, selection) %>%
            summarise(method="greedy", Key=pow.str, Value=mean(pow_greedy>0), SE=sd(pow_greedy>0)/sqrt(n()), N=n())
        df.pow <- df.pow %>%
            rbind(df.pow.greedy) %>%
            filter(method %in% method.values.tmp) %>%
            mutate(Method = factor(method, method.values.tmp, method.labels.tmp), Target="Closed testing")
    } else {
        df.pow <- df.pow %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels), Target="Closed testing")
    }
    df <- rbind(df.pow, df.lb) %>%
        filter(Key %in% key.values) %>%
        mutate(Key = factor(Key, key.values, key.labels))
    df.fdr.lb <- results %>%
        filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier==plot.classifier) %>%
        pivot_longer(c("disc_bh"), names_to="method", values_to="discoveries") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method, tune_size, selection) %>%
        summarise(Median=median(discoveries), Quantile.90=quantile(discoveries, 0.9), SE=sd(discoveries)/sqrt(n()), N=n()) %>%
        pivot_longer(c("Median", "Quantile.90"), names_to="Key", values_to="Value") %>%
        mutate(method = "lb_auto", Target="FDR") %>%
        filter(method %in% method.values.tmp) %>%
        mutate(Method = factor(method, method.values.tmp, method.labels.tmp))
    df.fdr.pow <- results %>%
        filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier==plot.classifier) %>%
        pivot_longer(c("disc_bh"), names_to="method", values_to="discoveries") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method, tune_size, selection) %>%
        summarise(Key=pow.str, Value=mean(discoveries>0), SE=sd(discoveries>0)/sqrt(n()), N=n()) %>%
        mutate(method = "lb_auto", Target="FDR") %>%
        filter(method %in% method.values) %>%
        mutate(Method = factor(method, method.values.tmp, method.labels.tmp))
    df.fdr <- rbind(df.fdr.pow, df.fdr.lb) %>%
        filter(Key %in% key.values) %>%
        mutate(Key = factor(Key, key.values, key.labels))
    df.range <- tibble(Key=c(pow.str, pow.str, "Median", "Median", "Quantile.90", "Quantile.90"), Value=c(0,1,0,n_test*max_out,0,n_test*max_out),
                       n_out=c(0,n_test*max_out,0,n_test*max_out,0,n_test*max_out), method="lb_auto") %>%
        mutate(Method = factor(method, method.values.tmp, method.labels.tmp), Target="Closed testing") %>%
        filter(Key %in% key.values) %>%
        mutate(Key = factor(Key, key.values, key.labels))
    df.ref <- tibble(Key=c(pow.str, pow.str, "Median", "Median", "Quantile.90", "Quantile.90"),
                     Value=c(0.1,0.1,0,n_test*max_out,0,n_test*max_out), n_out=c(0,n_test*max_out,0,n_test*max_out,0,n_test*max_out), method="lb_auto", Target="Closed testing") %>%
        mutate(Method = factor(method, method.values.tmp, method.labels.tmp)) %>%
        filter(Key %in% key.values) %>%
        mutate(Key = factor(Key, key.values, key.labels))
    if(include.BH) {
        pp <- df %>%
            ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method, linetype=Target)) +
            geom_point() +
            geom_line() +
            geom_point(data=df.range, aes(x=n_out, y=Value, color=Method, shape=Method), alpha=0) +
            geom_line(data=df.ref, aes(x=n_out, y=Value), linetype=2, color="black", alpha=1) +
            geom_line(data=df.fdr, aes(x=n_out, y=Value, linetype=Target), color="black", alpha=1) +
            facet_wrap(.~Key, labeller="label_value", scales="free") +
            scale_color_manual(values=color.scale.tmp) +
            scale_shape_manual(values=shape.scale.tmp) +
            scale_alpha_manual(values=alpha.scale.tmp) +
            scale_linetype_manual(values=c(1,3)) +
            xlab("True number of outliers") +
            ylab("") +
            theme_bw() +
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))
    } else {
        pp <- df %>%
            ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method, linetype=Target)) +
            geom_point() +
            geom_line() +
            geom_point(data=df.range, aes(x=n_out, y=Value, color=Method, shape=Method), alpha=0) +
            geom_line(data=df.ref, aes(x=n_out, y=Value), linetype=2, color="black", alpha=1) +
            facet_wrap(.~Key, labeller="label_value", scales="free") +
            scale_color_manual(values=color.scale.tmp) +
            scale_shape_manual(values=shape.scale.tmp) +
            scale_alpha_manual(values=alpha.scale.tmp) +
            scale_linetype_manual(values=c(1,3)) +
            xlab("True number of outliers") +
            ylab("") +
            theme_bw() +
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))
    }
    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_lhco_nt%d_lb_BH_%s_simple%s.pdf", fig.name, plot.n_train, include.BH, plot.simple)
        ggsave(pp, file=plot.file.1, height=2.25, width=8, units="in")
    } else {
        print(pp)
    }
    ## Companion table (only for the full / non-simple variant)
    if(!plot.simple) {
        df.table <- df %>%
            filter(Target=="Closed testing") %>%
            ungroup() %>%
            select(n_out, Key, Value, SE, Method) %>%
            group_by(n_out, Key) %>%
            mutate(Value_max = max(Value[Method!="Naive"])) %>%
            ungroup() %>%
            mutate(`Outliers`=n_out,
                   Value = ifelse(Key==pow.str, sprintf("%.2f (%.2f)", Value, SE),
                                  sprintf("%4d (%d)", round(Value), round(SE)))) %>%
            select(`Outliers`, Key, Value, Method) %>%
            pivot_wider(names_from = c(Method), values_from = Value, names_sort=TRUE) %>%
            arrange(Key, `Outliers`)
        tab <- df.table %>%
            select(-Key) %>%
            kable(format="latex", booktabs=TRUE, align = 'c', escape=FALSE) %>%
            pack_rows(index = c("90% Lower bound (median)" = length(table(df.table$`Outliers`)),
                                "90% Lower bound (90-th quantile)" = length(table(df.table$`Outliers`)),
                                "Power (global null)" = length(table(df.table$`Outliers`))
                                )) %>%
            ## Header span = number of method columns (everything except "Outliers").
            ## Auto-derived so it tracks whatever methods are present.
            add_header_above(c(" ", "Local testing procedure" = ncol(df.table) - 2))
        writeLines(tab, sprintf("tables/%s_lhco_nt%d_nc%d.tex", tab.name, plot.n_train, plot.n_cal))
    }
}

## ---- Run --------------------------------------------------------------------
results <- load_data("fig1")

n_test <- 10000
max_out <- 0.15

## Figure 1
make_plot_lower_bound_lhco(fig.name="fig1", plot.n_train=10000, plot.n_cal=2000, n_test=n_test, max_out=max_out,
                           include.BH=TRUE, plot.classifier="auto",
                           plot.simple=TRUE, plot.naive=TRUE, save.plot=TRUE)

## Figure A13 and Table A4 (Appendix): full version with all local tests
make_plot_lower_bound_lhco(fig.name="figA13", tab.name="tabA4", plot.n_train=10000, plot.n_cal=2000, n_test=n_test, max_out=max_out,
                           include.BH=TRUE, plot.classifier="auto",
                           plot.simple=FALSE, plot.naive=TRUE, save.plot=TRUE)
