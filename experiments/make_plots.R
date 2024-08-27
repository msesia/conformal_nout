options(width = 300)

library(tidyverse)
library(latex2exp)
library(RColorBrewer)
library(kableExtra)

plot.synthetic.0 <- TRUE
plot.synthetic.1_2 <- TRUE
plot.synthetic.1_2_selection <- TRUE
plot.data.selection <- TRUE
plot.synthetic.3 <- TRUE
plot.data.setup.4 <- TRUE
plot.data.2 <- TRUE
plot.lehmann <- FALSE
plot.lehmann.new <- FALSE
plot.lhco_1 <- FALSE
plot.lhco_selection <- TRUE

load_data <- function(setup) {
    idir <- sprintf("results_hpc/setup%d", setup)
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}

init_settings <- function(idx.exclude=NULL, names_ACODE=FALSE) {
    cbPalette <<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#6e57d2", "red")
    method.values <<- c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_lmp", "lb_auto") # "lb_wmw_k3", "lb_wmw_k4",
    if(names_ACODE) {
        method.labels <<- c("ACODE (Simes)", "ACODE (Storey-Simes)", "ACODE (Fisher)", "ACODE (WMW)", "ACODE (LMP, Lehmann k=3)", "ACODE (LMP,G-hat)", "ACODE (Adaptive)")
    } else {
        method.labels <<- c("Simes", "Storey-Simes", "Fisher", "WMW", "LMPI (Lehmann, k=3)", "LMPI (G-hat)", "Adaptive") #"WMW (k=3)",
    }
    classifier.values <<- c("occ-auto", "bc-auto", "auto")
    classifier.labels <<- c("One-Class", "Binary", "Automatic")
    color.scale <<- cbPalette[c(1,1,4,3,6,7,8)]
    shape.scale <<- c(2,6,3,1,0,9,8)
    alpha.scale <<- c(0.75,0.75,0.75,0.75,0.75,0.75,1)
    data.values <<- c("pendigits", "creditcard", "cover", "shuttle", "mammography", "aloi")
    data.labels <<- c("Pendigits", "Creditcard", "Covertype", "Shuttle", "Mammography", "ALOI")
    if(length(idx.exclude)>0) {
        ## Exclude Storey-simes
        method.values <<- method.values[-idx.exclude]
        method.labels <<- method.labels[-idx.exclude]
        color.scale <<- color.scale[-idx.exclude]
        shape.scale <<- shape.scale[-idx.exclude]
        alpha.scale <<- alpha.scale[-idx.exclude]
    }
}

if(plot.synthetic.0) {

    make_plot_lower_bound_proportion <- function(setup, reload=FALSE, plot.quantile=0.5) {
        init_settings(idx.exclude=NULL)
        if(reload) {
            results <- load_data(setup)
        }
        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_lmp"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(LB=quantile(lower_bound, plot.quantile), LB.se=sd(lower_bound)/sqrt(n()))
        df <- summary %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
        x.max <- max(df$n_out)
        pp <- df %>%
            ggplot(aes(x=n_out, y=LB, color=Method, shape=Method, alpha=Method)) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se)), width=0.01, alpha=0.5) +
            geom_abline(slope=1, intercept=0, linetype=2) +
            facet_grid(.~Classifier, labeller="label_both") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("True number of outliers") +
            ylab("90% lower bound") +
            xlim(0,x.max) +
            ylim(0,NA) +
            theme_bw() +
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))
        plot.file.1 <- sprintf("figures/synthetic%d_lower_bound_q%s.pdf", setup, plot.quantile)
        ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
    }

    make_plot_power_proportion <- function(setup, reload=FALSE, plot.quantile=0.5, plot.alpha=0.1) {
        init_settings(idx.exclude=NULL)
        if(reload) {
            results <- load_data(setup)
        }
        summary <- results %>%
            pivot_longer(c("pval_simes", "pval_storey_simes", "pval_fisher", "pval_auto", "pval_wmw_k2", "pval_wmw_k3", "pval_lmp"),
                         names_to="method", values_to="pval") %>%
            mutate(method = str_replace(method, "^pval_", "lb_")) %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(Power=mean(pval<=plot.alpha), Power.se=sd(pval<=plot.alpha)/sqrt(n()), N=n())
        df <- summary %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
        x.max <- max(df$n_out)
        pp <- df %>%
            ggplot(aes(x=n_out, y=Power, color=Method, shape=Method, alpha=Method)) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin=(Power-2*Power.se), ymax=(Power+2*Power.se)), width=0.01, alpha=0.5) +
            geom_abline(slope=0, intercept=plot.alpha, linetype=2) +
            facet_grid(.~Classifier, labeller="label_both") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("True number of outliers") +
            ylab("90% lower bound") +
            xlim(0,x.max) +
            ylim(0,NA) +
            theme_bw()
        pp
#        plot.file.1 <- sprintf("figures/synthetic%d_lower_bound_q%s.pdf", setup, plot.quantile)
#        ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
    }

    make_combined_plot_proportion <- function(setup, reload=FALSE, plot.quantile=0.5, plot.alpha=0.1) {
        init_settings(idx.exclude=NULL)
        if(reload) {
            results <- load_data(setup)
        }
        ## Lower Bound Summary
        lb_summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_lmp"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(Value=quantile(lower_bound, plot.quantile), Value.se=sd(lower_bound)/sqrt(n()), .groups = 'drop') %>%
            mutate(metric = "Lower Bound")
        ## Power Summary
        power_summary <- results %>%
            pivot_longer(c("pval_simes", "pval_storey_simes", "pval_fisher", "pval_auto", "pval_wmw_k2", "pval_wmw_k3", "pval_lmp"),
                         names_to="method", values_to="pval") %>%
            mutate(method = str_replace(method, "^pval_", "lb_")) %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(Value=mean(pval<=plot.alpha), Value.se=sd(pval<=plot.alpha)/sqrt(n()), .groups = 'drop') %>%
            mutate(metric = "Power (global null)")
        ## Combine summaries
        combined_summary <- bind_rows(lb_summary, power_summary)
        df <- combined_summary %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
        x.max <- max(df$n_out)
        ## Plot
        data.ref.pow <- tibble(n_out=c(0, max(df$n_out)), metric="Power (global null)", Value=c(0.1,0.1), Method="Adaptive") %>%
            mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound")))            
        data.ref.lb <- tibble(n_out=c(0, max(df$n_out)), metric="Lower Bound", Value=c(0, max(df$n_out)), Method="Adaptive")  %>%
        mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound")))            
        data.lim.pow <- tibble(n_out=c(0, 0), metric="Power (global null)", Value=c(0,1), Method="Adaptive")  %>%
        mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound")))            
        data.lim.lb <- tibble(n_out=c(0, 0), metric="Lower Bound", Value=c(0, max(df$n_out)), Method="Adaptive")  %>%
        mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound")))            
        pp <- df %>%
            mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound"))) %>%
            ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method)) +
            geom_point() +
            geom_line() +
##            geom_errorbar(aes(ymin=(Value-2*Value.se), ymax=(Value+2*Value.se)), width=0.01, alpha=0.5) +
            geom_line(data = data.ref.pow, linetype=2, color="black", show_guide = FALSE) +
            geom_line(data = data.ref.lb, linetype=2, color="black", show_guide = FALSE) +
            geom_point(data = data.lim.pow, alpha=0) +
            geom_point(data = data.lim.lb, alpha=0) +
            facet_grid(metric ~ Classifier, scales = "free_y", labeller = labeller(Classifier = label_both, metric = label_value)) +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("True number of outliers") +
            ylab("") +
            xlim(0,x.max) +
#            ylim(0,NA) +
            theme_bw() +
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))
        plot.file.1 <- sprintf("figures/synthetic%d_combined_plot_q%s_alpha%s.pdf", setup, plot.quantile, plot.alpha)
        ggsave(pp, file=plot.file.1, height=3.5, width=7, units="in")
    }

    results <- load_data(0)
   
    make_plot_lower_bound_proportion(0, plot.quantile=0.9, reload=TRUE)
    ##make_plot_power_proportion(0, plot.quantile=0.9, reload=TRUE)

    make_combined_plot_proportion(0, plot.quantile=0.9, reload=FALSE)

}


if(plot.synthetic.1_2) {

    make_plot_lower_bound_proportion <- function(setup, reload=FALSE, plot.quantile=0.5) {
        init_settings(idx.exclude=c(6))
        if(reload) {
            results <- load_data(setup)
        }
        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(LB=quantile(lower_bound, plot.quantile), LB.se=sd(lower_bound)/sqrt(n()))
        df <- summary %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
        x.max <- max(df$n_out)
        pp <- df %>%
            ggplot(aes(x=n_out, y=LB, color=Method, shape=Method, alpha=Method)) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se)), width=0.01, alpha=0.5) +
            geom_abline(slope=1, intercept=0, linetype=2) +
            facet_grid(.~Classifier, labeller="label_both") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("True number of outliers") +
            ylab("90% lower bound") +
            xlim(0,x.max) +
            ylim(0,x.max) +
            theme_bw() +
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))            
        plot.file.1 <- sprintf("figures/synthetic%d_lower_bound_q%s.pdf", setup, plot.quantile)
        ggsave(pp, file=plot.file.1, height=2.25, width=7.5, units="in")
    }

    make_plot_lower_bound_proportion(1, plot.quantile=0.5, reload=TRUE)
    make_plot_lower_bound_proportion(2, plot.quantile=0.5, reload=TRUE)

    make_plot_lower_bound_proportion(1, plot.quantile=0.9, reload=TRUE)
    make_plot_lower_bound_proportion(2, plot.quantile=0.9, reload=TRUE)

}

if(plot.synthetic.1_2_selection) {

    make_plot_lower_bound_proportion_sel <- function(setup, reload=FALSE) {
        init_settings(idx.exclude=c(6))
        if(reload) {
            results <- load_data(setup)
        }
        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier, method) %>%
            summarise(LB=median(lower_bound), LB.se=sd(lower_bound)/sqrt(n()))
        df <- summary %>%
            filter(prop_out %in% c(0,0.2,0.5)) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test))
        df.ref <- results %>%
            filter(prop_out %in% c(0,0.2,0.5)) %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier) %>%
            summarise(LB=median(n_out_sel)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test))
        df.ghost <- tibble(prop_out = c(0,0), selected_num=c(0,0), n_test=c(1000,1000), LB=c(0,1)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test))
        x.min <- min(df$selected_num/df$n_test)
        x.max <- max(df$selected_num/df$n_test)
        pp <- df %>%
            ggplot(aes(x=selected_num/n_test, y=LB)) +
            geom_point(aes(color=Method, shape=Method, alpha=Method)) +
            geom_line(aes(color=Method, alpha=Method)) +
            geom_line(data=df.ref, linetype=2, aes(x=selected_num/n_test, y=LB), color="black", alpha=1) +
            geom_point(data=df.ghost, aes(x=selected_num/n_test, y=LB), alpha=0) +
##            geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se), color=Method), width=0.01, alpha=0.5) +
            facet_wrap(.~N_out, labeller="label_value", scale="free", nrow=1) +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("Proportion of top scores in selected test set") +
            ylab("90% lower bound") +
    scale_x_continuous(trans='log10', limits=c(0.01,1)) +
            theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))                
        plot.file.1 <- sprintf("figures/synthetic%d_lower_bound_sel.pdf", setup)
        ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
    }

    make_plot_lower_bound_proportion_sel(1001, reload=TRUE)

    make_plot_lower_bound_proportion_sel(1002, reload=TRUE)

}


if(plot.synthetic.3) {

    make_plot_lower_bound_mixture <- function(plot.p, plot.n_train, plot.signal, reload=FALSE) {
        init_settings(idx.exclude=c(6))
        if(reload) {
            results <- load_data(3)
        }
        summary <- results %>%
            separate(Data, into = c("Data","Mixture"), sep = "-") %>%
            mutate(Mixture=parse_number(Mixture)) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"), names_to="method", values_to="lower_bound") %>%
            group_by(Data, Mixture, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(LB=mean(lower_bound), LB.se=sd(lower_bound)/sqrt(n()))
        df <- summary %>%
            filter(p==plot.p, n_train==plot.n_train, Signal==plot.signal) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
        pp <- df %>%
            ggplot(aes(x=Mixture, y=LB/n_test, color=Method, shape=Method, alpha=Method)) +
            geom_point() +
            geom_line() +
#            geom_errorbar(aes(ymin=(LB-2*LB.se)/n_test, ymax=(LB+2*LB.se)/n_test), width=0.01, alpha=0.5) +
            facet_grid(.~Classifier, labeller="label_both") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("Proportion of data points from binomial distribution") +
            ylab("90% lower bound") +
#            xlim(NA,0.5) +
#            ylim(NA,0.5) +
            theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))                
        plot.file.1 <- sprintf("figures/synthetic3_p%d_n%d_s%.2f_lower_bound.pdf", plot.p, plot.n_train, plot.signal)
        ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
    }

    make_plot_lower_bound_mixture(plot.p=100, plot.n_train=1000, plot.signal=3, reload=TRUE)

}

if(plot.data.setup.4) {

    make_plot_data_4 <- function(plot.data, plot.n_train, plot.n_cal, plot.n_test, reload=FALSE, save.plot=True) {
        init_settings(idx.exclude=c(6))
        if(reload){
            results <- load_data(4)
        }
        key.values <- c("LB.50", "LB.90", "Power")
        key.labels <- c("LB (median)", "LB (90th q.)", "Power (global)")
        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(Value_LB.50=median(lower_bound), Value_LB.90=quantile(lower_bound, 0.9), Value_Power=mean(lower_bound>0),
                      SE_LB.50=sd(lower_bound)/sqrt(n()), SE_LB.90=sd(lower_bound)/sqrt(n()), SE_Power=sd(lower_bound>0)/sqrt(n()), N=n()) %>%
            pivot_longer(c("Value_LB.50", "Value_LB.90", "Value_Power", "SE_LB.50", "SE_LB.90", "SE_Power"),
                         names_to = c(".value", "Key"), names_sep = "_")
        df <- summary %>%
            filter(Data==plot.data, n_train==plot.n_train, n_cal==plot.n_cal, n_test==plot.n_test) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels)) %>%
            mutate(Key = factor(Key, key.values, key.labels))
        x.max <- max(df$n_out)
        df.ref <- tibble(Key=c("LB.50", "LB.50", "LB.90", "LB.90", "Power", "Power"),
                         n_out=c(0,x.max,0,x.max,0,x.max), Value=c(0,x.max,0,x.max,0.1,0.1),
                         method="lb_auto", Target="Closed testing") %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Key = factor(Key, key.values, key.labels))
        df.range <- tibble(Key=c("LB.50", "LB.50", "LB.90", "LB.90", "Power", "Power"),
                         n_out=c(0,x.max,0,x.max,0,x.max), Value=c(0,x.max,0,x.max,0,1),
                         method="lb_auto", Target="Closed testing") %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Key = factor(Key, key.values, key.labels))
        pp <- df %>%
            ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method)) +
            geom_point() +
            geom_line() +
            ##geom_errorbar(aes(ymin=(Value-2*SE), ymax=(Value+2*SE)), width=0.1, alpha=0.5) +
            geom_line(data=df.ref, aes(x=n_out, y=Value), linetype=2, color="black", alpha=1) +
            geom_point(data=df.range, aes(x=n_out, y=Value), color="black", alpha=0) +
            facet_grid(Key~Classifier, labeller=labeller(Classifier = label_both, Key = label_value), scales="free") +
            xlab("True number of outliers") +
            ylab("") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlim(0,x.max) +
            theme_bw() +
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))
        pp        
        if(save.plot) {
            plot.file.1 <- sprintf("figures/setup4_%s_%d_%d_%d.pdf", plot.data, plot.n_train, plot.n_cal, plot.n_test)
            ggsave(pp, file=plot.file.1, height=3.5, width=7, units="in")
        } else {
            print(pp)
        }
    }


    results <- load_data(4)

    data.list <- c("creditcard", "pendigits", "cover", "shuttle", "mammography", "aloi")
    ##data.list <- c("creditcard")

    plot.n_train <- 1000
    plot.n_cal <- 200
    for(plot.n_test in c(100)) {
        for(plot.data in data.list) { #
            make_plot_data_4(plot.data, plot.n_train, plot.n_cal, plot.n_test, reload=FALSE, save.plot=TRUE)
        }
    }

}


if(plot.data.1) {

    make_plot_lower_bound_data_5 <- function(plot.data, plot.n_train, plot.classifier, reload=FALSE, save.plot=True) {
        if(reload){
            results <- load_data(5)
        }
        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(LB=mean(lower_bound), LB.se=sd(lower_bound)/sqrt(n()), N=n())
        df <- summary %>%
            filter(Data==plot.data, n_train==plot.n_train, Classifier==plot.classifier) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
        pp <- df %>%
            ggplot(aes(x=n_test, y=LB/n_test, color=Method, shape=Method, alpha=Method)) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin=(LB-2*LB.se)/n_test, ymax=(LB+2*LB.se)/n_test), width=0.1, alpha=0.5) +
            geom_abline(slope=1, intercept=0, linetype=2) +
            facet_grid(prop_out~n_cal, labeller="label_both") +
            xlab("Size of test set") +
            ylab("90% lower bound") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            theme_bw()
        if(save.plot) {
            plot.file.1 <- sprintf("figures/setup5_%s_lb_%d_%s.png", plot.data, plot.n_train, plot.classifier)
            ggsave(pp, file=plot.file.1, height=5, width=7, units="in")
        } else {
            print(pp)
        }
    }

    make_plot_power_data_5 <- function(plot.data, plot.n_train, plot.classifier, reload=FALSE, save.plot=True) {
        if(reload){
            results <- load_data(5)
        }
        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3"),
                         names_to="method", values_to="lower_bound") %>%
            mutate(Reject=(lower_bound>0)) %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(Power=mean(Reject), Power.se=sd(Reject)/sqrt(n()), N=n())
        df <- summary %>%
            filter(Data==plot.data, n_train==plot.n_train, Classifier==plot.classifier) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
        x.max <- max(df$n_out)
        pp <- df %>%
            ggplot(aes(x=n_test, y=Power, color=Method, shape=Method, alpha=Method)) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin=(Power-2*Power.se), ymax=(Power+2*Power.se)), width=0.1, alpha=0.5) +
            geom_abline(slope=0, intercept=0.1, linetype=2) +
            facet_grid(prop_out~n_cal, labeller="label_both") +
            xlab("Size of test set") +
            ylab("Power") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            ylim(0, NA) +
            theme_bw()
        if(save.plot) {
            plot.file.1 <- sprintf("figures/setup5_%s_power_%d_%s.png", plot.data, plot.n_train, plot.classifier)
            ggsave(pp, file=plot.file.1, height=5, width=7, units="in")
        } else {
            print(pp)
        }
    }

    results <- load_data(5)

    data.list <- c("creditcard", "pendigits", "cover")
    #data.list <- c("creditcard")

    plot.n_train <- 1000
    for(plot.data in data.list) { #
        for(plot.classifier in c("occ-if", "auto")) { #
            make_plot_lower_bound_data_5(plot.data, plot.n_train, plot.classifier, reload=FALSE, save.plot=TRUE)
            make_plot_power_data_5(plot.data, plot.n_train, plot.classifier, reload=FALSE, save.plot=TRUE)
        }
    }

 }




if(plot.data.selection) {

    make_plot_lower_bound_data_sel <- function(setup, plot.quantile=0.5, reload=FALSE) {
        init_settings(idx.exclude=c(6))
        if(reload) {
            results <- load_data(setup)
        }
        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3"), names_to="method", values_to="lower_bound") %>% # , "lb_wmw_k4"
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier, method) %>%
            summarise(LB=quantile(lower_bound, plot.quantile), LB.se=sd(lower_bound)/sqrt(n()))
        df <- summary %>%
            filter(prop_out %in% c(0,0.2,0.5)) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
            mutate(Data = factor(Data, data.values, data.labels))
        df.ref <- results %>%
            filter(prop_out %in% c(0,0.2,0.5)) %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier) %>%
            summarise(LB=median(n_out_sel)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
            mutate(Data = factor(Data, data.values, data.labels))
        df.ghost <- tibble(prop_out = c(0,0), selected_num=c(0,0), n_test=c(1000,1000), LB=c(0,1)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test))
        x.min <- min(df$selected_num/df$n_test)
        x.max <- max(df$selected_num/df$n_test)
        pp <- df %>%
            ggplot(aes(x=selected_num/n_test, y=LB)) +
            geom_point(aes(color=Method, shape=Method, alpha=Method)) +
            geom_line(aes(color=Method, alpha=Method)) +
            geom_line(data=df.ref, linetype=2, aes(x=selected_num/n_test, y=LB), color="black", alpha=1) +
            geom_point(data=df.ghost, aes(x=selected_num/n_test, y=LB), alpha=0) +
                                        #            geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se), color=Method), width=0.01, alpha=0.5) +
            facet_wrap(.~Data, labeller="label_value", scale="free", nrow=2) +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("Proportion of top scores in selected test set") +
            ylab("90% lower bound") +
            scale_x_continuous(trans='log10', limits=c(0.01,1)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        plot.file.1 <- sprintf("figures/data%d_lower_bound_sel_q%s.pdf", setup, plot.quantile)
        ggsave(pp, file=plot.file.1, height=3.5, width=7, units="in")
    }


    results <- load_data(1005)
    make_plot_lower_bound_data_sel(1005, plot.quantile=0.5, reload=FALSE)
    make_plot_lower_bound_data_sel(1005, plot.quantile=0.9, reload=FALSE)

}



if(FALSE) {



    make_plot_lower_bound_1 <- function(plot.p, plot.n_train, plot.signal, reload=FALSE) {
        if(reload) {
            results <- load_data(1)
        }

        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"), names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(LB=mean(lower_bound), LB.se=sd(lower_bound)/sqrt(n()))

        cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#        method.values = c("lb_bh", "lb_simes", "lb_storey_simes", "lb_wmw", "lb_wmw_k3")
#        method.labels = c("d BH", "d Simes", "d ASimes", "d WMW", "d WMW (k3)")
        method.values = c("lb_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_auto")
        method.labels = c("d Simes", "d Fisher", "d WMW (k2)", "d WMW (k3)", "d Auto")

        classifier.values = c("occ", "bc", "auto")
        classifier.labels = c("One-Class", "Binary", "Auto")

#        color.scale <- cbPalette[c(1,3,6,7,8)]
#        color.scale <- cbPalette[c(1,3,4,7,7)]
#        shape.scale <- c(5,1,16,0,15,2)

        df <- summary %>%
            filter(p==plot.p, n_train==plot.n_train, Signal==plot.signal) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))

        pp <- df %>%
            ggplot(aes(x=n_out/n_test, y=LB/n_test, color=Method, shape=Method)) +
            geom_point() +
            geom_line() +
#            geom_errorbar(aes(ymin=(LB-2*LB.se)/n_test, ymax=(LB+2*LB.se)/n_test), width=0.01, alpha=0.5) +
            geom_abline(slope=1, intercept=0, linetype=2) +
            facet_grid(.~Classifier, labeller="label_both") +
#            scale_color_manual(values=color.scale) +
#            scale_shape_manual(values=shape.scale) +
            xlab("True proportion of outliers") +
            ylab("90% lower bound") +
            xlim(NA,0.5) +
            ylim(NA,0.5) +
            theme_bw()
        plot.file.1 <- sprintf("figures/synthetic1_p%d_n%d_s%.2f_lower_bound.png", plot.p, plot.n_train, plot.signal)
        ggsave(pp, file=plot.file.1, height=2.5, width=7, units="in")
    }

    make_plot_lower_bound_1(plot.p=1000, plot.n_train=1000, plot.signal=0.7, reload=TRUE)



    make_plot_lower_bound_2 <- function(plot.p, plot.n_train, plot.signal, reload=FALSE) {
        if(reload) {
            results <- load_data(2)
        }

        summary <- results %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"), names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(LB=mean(lower_bound), LB.se=sd(lower_bound)/sqrt(n()))

        cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#        method.values = c("lb_bh", "lb_simes", "lb_storey_simes", "lb_wmw", "lb_wmw_k3")
#        method.labels = c("d BH", "d Simes", "d ASimes", "d WMW", "d WMW (k3)")
        method.values = c("lb_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_auto")
        method.labels = c("d Simes", "d Fisher", "d WMW (k2)", "d WMW (k3)", "d Auto")

        classifier.values = c("occ", "bc", "auto")
        classifier.labels = c("One-Class", "Binary", "Auto")

#        color.scale <- cbPalette[c(1,3,6,7,8)]
#        color.scale <- cbPalette[c(1,3,4,7,7)]
#        shape.scale <- c(5,1,16,0,15,2)

        df <- summary %>%
            filter(p==plot.p, n_train==plot.n_train, Signal==plot.signal) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))

        pp <- df %>%
            ggplot(aes(x=n_out/n_test, y=LB/n_test, color=Method, shape=Method)) +
            geom_point() +
            geom_line() +
#            geom_errorbar(aes(ymin=(LB-2*LB.se)/n_test, ymax=(LB+2*LB.se)/n_test), width=0.01, alpha=0.5) +
            geom_abline(slope=1, intercept=0, linetype=2) +
            facet_grid(.~Classifier, labeller="label_both") +
#            scale_color_manual(values=color.scale) +
#            scale_shape_manual(values=shape.scale) +
            xlab("True proportion of outliers") +
            ylab("90% lower bound") +
            xlim(NA,0.5) +
            ylim(NA,0.5) +
            theme_bw()
        plot.file.1 <- sprintf("figures/synthetic2_p%d_n%d_s%d_lower_bound.png", plot.p, plot.n_train, plot.signal)
        ggsave(pp, file=plot.file.1, height=2.5, width=7, units="in")
    }

    make_plot_lower_bound_2(plot.p=10, plot.n_train=1000, plot.signal=6, reload=TRUE)


    make_plot_lower_bound_3 <- function(plot.p, plot.n_train, plot.signal, reload=FALSE) {
        if(reload) {
            results <- load_data(3)
        }

        summary <- results %>%
            separate(Data, into = c("Data","Mixture"), sep = "-") %>%
            mutate(Mixture=parse_number(Mixture)) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"), names_to="method", values_to="lower_bound") %>%
            group_by(Data, Mixture, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
            summarise(LB=mean(lower_bound), LB.se=sd(lower_bound)/sqrt(n()))

        cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#        method.values = c("lb_bh", "lb_simes", "lb_storey_simes", "lb_wmw", "lb_wmw_k3")
#        method.labels = c("d BH", "d Simes", "d ASimes", "d WMW", "d WMW (k3)")
        method.values = c("lb_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_auto")
        method.labels = c("d Simes", "d Fisher", "d WMW (k2)", "d WMW (k3)", "d Auto")

        classifier.values = c("occ", "bc", "auto")
        classifier.labels = c("One-Class", "Binary", "Auto")

#        color.scale <- cbPalette[c(1,3,6,7,8)]
#        color.scale <- cbPalette[c(1,3,4,7,7)]
#        shape.scale <- c(5,1,16,0,15,2)

        df <- summary %>%
            filter(p==plot.p, n_train==plot.n_train, Signal==plot.signal) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))

        pp <- df %>%
            ggplot(aes(x=Mixture, y=LB/n_test, color=Method, shape=Method)) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin=(LB-2*LB.se)/n_test, ymax=(LB+2*LB.se)/n_test), width=0.01, alpha=0.5) +
            geom_abline(slope=0, intercept=0.5, linetype=2) +
            facet_grid(.~Classifier, labeller="label_both") +
#            scale_color_manual(values=color.scale) +
#            scale_shape_manual(values=shape.scale) +
            xlab("Mixture") +
            ylab("90% lower bound") +
#            xlim(NA,0.5) +
#            ylim(NA,0.5) +
            theme_bw()
        plot.file.1 <- sprintf("figures/synthetic3_p%d_n%d_s%.2f_lower_bound.png", plot.p, plot.n_train, plot.signal)
        ggsave(pp, file=plot.file.1, height=2.5, width=7, units="in")
    }

    make_plot_lower_bound_3(plot.p=100, plot.n_train=1000, plot.signal=3, reload=TRUE)


    make_plot_discoveries <- function(plot.model) {

        summary <- results %>%
            pivot_longer(
                c("fdp_bh", "fdp_sbh", "power_bh", "power_sbh"),
                names_to = c(".value", "method"),
                names_pattern = "(.*)_(.*)"
            ) %>%
            mutate(FDR=fdp, Power=power) %>%
            pivot_longer(c("FDR", "Power"), names_to="metric", values_to="value") %>%
            group_by(alpha, prop_out, n_test, model_name, n_ref, n_cal, n_out, method, metric) %>%
            summarise(Value=mean(value), SE=sd(value)/sqrt(n()))

        cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "grey20")
        method.values = c("bh", "sbh")
        method.labels = c("BH", "ABH")
        method.values = c("bh")
        method.labels = c("BH")
        color.scale <- cbPalette[c(1,9)]
        shape.scale <- c(5,5)

        df.nominal <- tibble(metric="FDR", Value=0.1, prop_out=0)

        df <- summary %>%
            filter(model_name==plot.model) %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels))

        pp <- df %>%
            ggplot(aes(x=n_out, y=Value, color=Method, shape=Method)) +
            geom_point() +
            geom_line() +
#            geom_errorbar(aes(ymin=Value-2*SE, ymax=Value+2*SE), width=0.1, alpha=0.5) +
            facet_grid(n_ref~metric) +
            geom_hline(data=df.nominal, aes(yintercept=Value), linetype=2) +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            xlab("True number of signal events") +
            ylim(0,1) +
            theme_bw()

        plot.file.1 <- sprintf("figures/lhco_discoveries_%s.png", plot.model)
        ggsave(pp, file=plot.file.1, height=3, width=5, units="in")

    }

    make_plot_discoveries("bc-abc")


}


load_data <- function(setup) {
    idir <- sprintf("results_hpc/setup%d", setup)
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}


make_plot_power_data_1 <- function(plot.data, plot.n_train, plot.n_cal, plot.n_test, reload=FALSE, save.plot=True) {
    if(reload){
        results <- load_data(4)
    }

    summary <- results %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
                     names_to="method", values_to="lower_bound") %>%
        mutate(Reject=(lower_bound>0)) %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
        summarise(Power=mean(Reject), Power.se=sd(Reject)/sqrt(n()), N=n())

    cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        #        method.values = c("lb_bh", "lb_simes", "lb_storey_simes", "lb_wmw", "lb_wmw_k3")
                                        #        method.labels = c("d BH", "d Simes", "d ASimes", "d WMW", "d WMW (k3)")
    method.values = c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_auto")
    method.labels = c("d Simes", "d Simes (Storey)", "d Fisher", "d WMW (k2)", "d WMW (k3)", "d Auto")

    classifier.values = c("occ", "bc", "auto")
    classifier.labels = c("One-Class", "Binary", "Auto")

    df <- summary %>%
        filter(Data==plot.data, n_train==plot.n_train, n_cal==plot.n_cal, n_test==plot.n_test) %>%
        filter(method %in% method.values) %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))

    pp <- df %>%
        ggplot(aes(x=n_out, y=Power, color=Method, shape=Method)) +
        geom_point() +
        geom_line() +
        geom_errorbar(aes(ymin=(Power-2*Power.se), ymax=(Power+2*Power.se)), width=0.1, alpha=0.5) +
        geom_abline(slope=0, intercept=0.1, linetype=2) +
        facet_grid(.~Classifier, labeller="label_both") +
                                        #            scale_color_manual(values=color.scale) +
                                        #            scale_shape_manual(values=shape.scale) +
        xlab("True number of outliers") +
        ylab("Power") +
                                        #            xlim(NA,0.5) +
                                        #            ylim(NA,0.5) +
        theme_bw()

    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_power_%d_%d_%d.png", plot.data, plot.n_train, plot.n_cal, plot.n_test)
        ggsave(pp, file=plot.file.1, height=3, width=5, units="in")
    } else {
        print(pp)
    }
}

make_plot_lower_bound_data_1 <- function(plot.data, plot.n_train, plot.n_cal, plot.n_test, reload=FALSE, save.plot=True) {
    if(reload){
        results <- load_data(4)
    }

    summary <- results %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
                     names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
        summarise(LB=mean(lower_bound), LB.se=sd(lower_bound)/sqrt(n()), N=n())

    cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        #        method.values = c("lb_bh", "lb_simes", "lb_storey_simes", "lb_wmw", "lb_wmw_k3")
                                        #        method.labels = c("d BH", "d Simes", "d ASimes", "d WMW", "d WMW (k3)")
    method.values = c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_auto")
    method.labels = c("d Simes", "d Simes (Storey)", "d Fisher", "d WMW (k2)", "d WMW (k3)", "d Auto")

    classifier.values = c("occ", "bc", "auto")
    classifier.labels = c("One-Class", "Binary", "Auto")

    df <- summary %>%
        filter(Data==plot.data, n_train==plot.n_train, n_cal==plot.n_cal, n_test==plot.n_test) %>%
        filter(method %in% method.values) %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))

    pp <- df %>%
        ggplot(aes(x=n_out, y=LB, color=Method, shape=Method)) +
        geom_point() +
        geom_line() +
        geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se)), width=0.1, alpha=0.5) +
        geom_abline(slope=1, intercept=0, linetype=2) +
        facet_grid(.~Classifier, labeller="label_both") +
                                        #            scale_color_manual(values=color.scale) +
                                        #            scale_shape_manual(values=shape.scale) +
        xlab("True number of outliers") +
        ylab("90% lower bound") +
                                        #            xlim(NA,0.5) +
                                        #            ylim(NA,0.5) +
        theme_bw()

    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_lb_%d_%d_%d.png", plot.data, plot.n_train, plot.n_cal, plot.n_test)
        ggsave(pp, file=plot.file.1, height=3, width=5, units="in")
    } else {
        print(pp)
    }
}


results <- load_data(4)

plot.n_train <- 500
plot.n_cal <- 500
for(plot.n_test in c(20, 100)) {
    for(plot.data in c("creditcard", "pendigits", "cover")) { #
        make_plot_power_data_1(plot.data, plot.n_train, plot.n_cal, plot.n_test, reload=FALSE, save.plot=TRUE)
        make_plot_lower_bound_data_1(plot.data, plot.n_train, plot.n_cal, plot.n_test, reload=FALSE, save.plot=TRUE)
    }
}

make_plot_lower_bound_data_1("pendigits", 199, 199, 20, reload=TRUE)


make_plot_lower_bound_data_test_size <- function(plot.data, plot.n_train, plot.n_cal, plot.prop_out=0.2, reload=FALSE, save.plot=True) {
    if(reload){
        results <- load_data(5)
    }

    summary <- results %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
                     names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
        summarise(LB=mean(lower_bound), LB.se=sd(lower_bound)/sqrt(n()), N=n())

    cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        #        method.values = c("lb_bh", "lb_simes", "lb_storey_simes", "lb_wmw", "lb_wmw_k3")
                                        #        method.labels = c("d BH", "d Simes", "d ASimes", "d WMW", "d WMW (k3)")
    method.values = c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_auto")
    method.labels = c("d Simes", "d Simes (Storey)", "d Fisher", "d WMW (k2)", "d WMW (k3)", "d Auto")

    classifier.values = c("occ", "bc", "auto")
    classifier.labels = c("One-Class", "Binary", "Auto")

    df <- summary %>%
        filter(Data==plot.data, n_train==plot.n_train, n_cal==plot.n_cal, prop_out==plot.prop_out) %>%
        filter(method %in% method.values) %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))

    pp <- df %>%
        ggplot(aes(x=n_out, y=LB, color=Method, shape=Method)) +
        geom_point() +
        geom_line() +
        geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se)), width=0.1, alpha=0.5) +
        geom_abline(slope=1, intercept=0, linetype=2) +
        facet_grid(.~Classifier, labeller="label_both") +
                                        #            scale_color_manual(values=color.scale) +
                                        #            scale_shape_manual(values=shape.scale) +
        xlab("True number of outliers") +
        ylab("90% lower bound") +
                                        #            xlim(NA,0.5) +
                                        #            ylim(NA,0.5) +
        theme_bw()

    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_lb_%d_%d_ntest_out%.2f.png", plot.data, plot.n_train, plot.n_cal, plot.prop_out)
        ggsave(pp, file=plot.file.1, height=3, width=5, units="in")
    } else {
        print(pp)
    }

}


results <- load_data(5)

plot.n_train <- 500
plot.n_cal <- 500
for(plot.prop_out in c(0.2)) {
    for(plot.data in c("pendigits", "creditcard", "cover")) {
        make_plot_lower_bound_data_test_size(plot.data, plot.n_train, plot.n_cal, plot.prop_out, reload=FALSE, save.plot=TRUE)
    }
}




if(plot.lhco_1) {

    make_plot_lower_bound_lhco <- function(reload=FALSE, plot.n_train, plot.n_cal, include.BH=FALSE,  plot.classifier="auto",
                                           plot.simple=FALSE, plot.naive=FALSE, save.plot=True, alpha=0.1) {
        if(plot.simple) {
            init_settings(idx.exclude=c(2,5,6))
        } else {
            init_settings(idx.exclude=c(6))
        }
        if(reload){
            results <- load_data(100)
        }
        pow.str <- "Power (global null)"
##        key.values <- c("Median", "Quantile.90", pow.str)
        key.values <- c("Median", pow.str)
##        key.labels <- c("Lower bound (median)", "Lower bound (90th quant.)", pow.str)
        key.labels <- c("Lower bound", pow.str)
        disc.str <- "Discoveries (10% FDR)"
        df.lb <- results %>%
            filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier==plot.classifier) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_auto"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method, tune_size, selection) %>%
            summarise(Median=median(lower_bound), Quantile.90=quantile(lower_bound, 0.9), SE=sd(lower_bound)/sqrt(n()), N=n()) %>%
            pivot_longer(c("Median", "Quantile.90"), names_to="Key", values_to="Value")
        if(plot.naive) {
            method.values.tmp <- c(method.values, "greedy")
            method.labels.tmp <- c(method.labels, "Cherry picking")
            color.scale.tmp <- c(color.scale, cbPalette[2])
            shape.scale.tmp <- c(shape.scale, 4)
            alpha.scale.tmp <- c(alpha.scale, 0.5)
            df.lb.greedy <- results %>%
                filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier!="auto") %>%
                pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_auto"),
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
            color.scale.tmp <- color.scale
            shape.scale.tmp <- shape.scale
            alpha.scale.tmp <- alpha.scale
            df.lb <- df.lb %>%
                filter(method %in% method.values) %>%
                mutate(Method = factor(method, method.values, method.labels), Target="Closed testing")
        }
        df.pow <- results %>%
            filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier==plot.classifier) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_auto"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method, tune_size, selection) %>%
            summarise(Key=pow.str, Value=mean(lower_bound>0), SE=sd(lower_bound>0)/sqrt(n()), N=n())
        if(plot.naive) {
            df.pow.greedy <- results %>%
                filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier!="auto") %>%
                pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_auto"),
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
        df.range <- tibble(Key=c(pow.str, pow.str, "Median", "Median", "Quantile.90", "Quantile.90"), Value=c(0,1,0,1500,0,1500),
                           n_out=c(0,1500,0,1500,0,1500), method="lb_auto") %>%
            mutate(Method = factor(method, method.values.tmp, method.labels.tmp), Target="Closed testing") %>%
            filter(Key %in% key.values) %>%
            mutate(Key = factor(Key, key.values, key.labels))
        df.ref <- tibble(Key=c(pow.str, pow.str, "Median", "Median", "Quantile.90", "Quantile.90"),
                         Value=c(0.1,0.1,0,1500,0,1500), n_out=c(0,1500,0,1500,0,1500), method="lb_auto", Target="Closed testing") %>%
            mutate(Method = factor(method, method.values.tmp, method.labels.tmp)) %>%
            filter(Key %in% key.values) %>%
            mutate(Key = factor(Key, key.values, key.labels))
        if(include.BH) {
            pp <- df %>%
                ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method, linetype=Target)) +
                geom_point() +
                geom_line() +
##                geom_errorbar(aes(ymin=(Value-2*SE), ymax=(Value+2*SE)), width=0.1, alpha=0.5) +
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
                ##geom_errorbar(aes(ymin=(Value-2*SE), ymax=(Value+2*SE)), width=0.1, alpha=0.5) +
                geom_point(data=df.range, aes(x=n_out, y=Value, color=Method, shape=Method), alpha=0) +
                geom_line(data=df.ref, aes(x=n_out, y=Value), linetype=2, color="black", alpha=1) +
                                        #                geom_line(data=df.fdr, aes(x=n_out, y=Value, linetype=Target), color="black", alpha=1) +
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
            plot.file.1 <- sprintf("figures/lhco_nt%d_lb_BH_%s_simple%s.pdf", plot.n_train, include.BH, plot.simple)
            ggsave(pp, file=plot.file.1, height=2.25, width=6, units="in")
        } else {
            print(pp)
        }
        ## Make table
        if(!plot.simple) {
            df.table <- df %>%
                filter(Target=="Closed testing") %>%
                ungroup() %>%
                select(n_out, Key, Value, SE, Method) %>%
                group_by(n_out, Key) %>%
                mutate(Value_max = max(Value[Method!="Naive"])) %>%
                                        #            mutate(Color = ifelse(Value==Value_max, "darkgreen", "black")) %>%
                ungroup() %>%
                mutate(`Outliers`=n_out,
                                        #                   Value = ifelse(Key==pow.str, sprintf("\\color{%s}{%.2f (%.2f)}", Color, Value, SE),
                                        #                                  sprintf("\\color{%s}{%4d (%d)}", Color, round(Value), round(SE)))) %>%
                       Value = ifelse(Key==pow.str, sprintf("%.2f (%.2f)", Value, SE),
                                      sprintf("%4d (%d)", round(Value), round(SE)))) %>%               
#                mutate(Key = factor(Key, key.values, key.labels)) %>%
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
                add_header_above(c(" ", "Local testing procedure" = 8))
            writeLines(tab, sprintf("tables/lhco_nt%d_nc%d.tex", plot.n_train, plot.n_cal))
        }
    }

    results <- load_data(100)


    make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=10000, plot.n_cal=2000, include.BH=TRUE, plot.classifier="auto", plot.simple=FALSE, plot.naive=TRUE, save.plot=TRUE)
    make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=10000, plot.n_cal=2000, include.BH=TRUE, plot.classifier="auto", plot.simple=TRUE, plot.naive=TRUE, save.plot=TRUE)

    make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=100000, plot.n_cal=2000, include.BH=TRUE, plot.classifier="bc-abc", plot.simple=FALSE, plot.naive=TRUE, save.plot=TRUE)


    #make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=100000, plot.n_cal=2000, include.BH=FALSE, save.plot=TRUE)
    #make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=10000, plot.n_cal=2000, include.BH=TRUE, save.plot=TRUE)
    #make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=100000, plot.n_cal=2000, include.BH=TRUE, save.plot=TRUE)



}

if(plot.lehmann.new) {

    load_data_la <- function(setup) {
        idir <- sprintf("results_hpc/lehmann%d", setup)
        ifile.list <- list.files(idir)
        results <- do.call("rbind", lapply(ifile.list, function(ifile) {
            df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
        }))
        return(results)
    }
    results <- load_data_la(91)

    plot.lehmann.new <- function(reload=FALSE) {
        if(reload){
            results <- load_data_la(91)
        }
        init_settings(idx.exclude=NULL)
        method.values <- c("Simes", "Storey.Simes", "Fisher", "WMW.k1", "WMW.k2", "WMW.k3")
        method.labels <- c("Simes", "Storey-Simes", "Fisher", "WMW (k=1)", "WMW (k=2)", "WMW (k=3)")
        dimensions.values <- c("m=100, n=100", "m=1000, n=100")
        summary <- results %>%
            group_by(Data, m, n, k, theta, Method) %>%
            summarise(Rejected=sum(Rejected), N=sum(N), Power=Rejected/N, Power.SE=sqrt(Power*(1-Power))/sqrt(N))
        plot.data <- "uniform"
        df <- summary %>%
            filter(Method %in% method.values) %>%
            mutate(Method = factor(Method, method.values, method.labels)) %>%
            mutate(Dimensions = sprintf("m=%d, n=%d", m, n), Dimensions=factor(Dimensions, dimensions.values)) %>%            
            mutate(Alterantive = sprintf("Lehmann's alternative, k=%d", k-1))
        pp <- df %>%
            ggplot(aes(x=theta, y=Power, color=Method, shape=Method)) +
            geom_point(alpha=0.75) +
            geom_line(alpha=0.75) +
                                        #        geom_errorbar(aes(ymin=Power-2*Power.SE, ymax=Power+2*Power.SE), alpha=0.5) +
            geom_hline(yintercept=0.1, linetype=2) +
            facet_grid(Dimensions~Alterantive, scales="free") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            ylim(0,NA) +
            xlim(NA,0.15) +
            xlab("Proportion of outliers") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        plot.file.1 <- sprintf("figures/lehmann_%s.pdf", plot.data)
        ggsave(pp, file=plot.file.1, height=4.5, width=7, units="in")
        ## Make table
        plot.m <- 100
        plot.n <- 100
        for(short.table in c(TRUE,FALSE)) {
            if(short.table) {
                theta.values <- c(0.10, 0.12, 0.14)
            } else {
                theta.values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15)
            }
            df.tab <- df %>%
                filter(m==plot.m, n==plot.n, theta %in% theta.values) %>%
                mutate(`Proportion of outliers`=theta) %>% select(-theta) %>%
                select(Dimensions, k, `Proportion of outliers`, Method, Power, Power.SE) %>%
                group_by(Dimensions, k, `Proportion of outliers`) %>%
                mutate(Max=max(Power)) %>%
                ungroup() %>%
                mutate(Power.str = sprintf("%.4f (%.4f)", Power, 2*Power.SE),
                       Power.str = ifelse(Power==Max, sprintf("\\textbf{%s}", Power.str), Power.str)) %>%
                select(k, `Proportion of outliers`, Method, Power.str) %>%       
                pivot_wider(names_from = c(Method), values_from = Power.str, names_sort=TRUE) %>%
                arrange(k, `Proportion of outliers`)
            tab <- df.tab %>%
                select(-k) %>%
                kable(format="latex", booktabs=TRUE, align = 'c', escape=FALSE) %>%                
                pack_rows("Lehmann's alternative with k=1", start_row = 1, end_row = cumsum(table(df.tab$k))[1]) %>%
                pack_rows("Lehmann's alternative with k=2", start_row = cumsum(table(df.tab$k))[1]+1, end_row = cumsum(table(df.tab$k))[2]) %>%
                pack_rows("Lehmann's alternative with k=3", start_row = cumsum(table(df.tab$k))[2]+1, end_row = cumsum(table(df.tab$k))[3]) %>%
                column_spec(c(1), width = "5em") %>%
                add_header_above(c(" "=1, "Testing procedure" = 6))                
            writeLines(tab, sprintf("tables/lehmann_short%s_m%d_n%d.tex", short.table, plot.m, plot.n))
        }
    }

    plot.lehmann.new(reload=FALSE)
            

    plot.lehmann.new.calibration <- function(reload=FALSE) {
        if(reload){
            results <- load_data_la(92)
        }
        init_settings(idx.exclude=NULL)
        method.values <- c("Simes", "Storey.Simes", "Fisher", "WMW.k1", "WMW.k2", "WMW.k3")
        method.labels <- c("Simes", "Storey-Simes", "Fisher", "WMW (k=1)", "WMW (k=2)", "WMW (k=3)")
        dimensions.values <- c("m=100, n=100", "m=1000, n=100")
        summary <- results %>%
            group_by(Data, m, n, k, theta, Method) %>%
            summarise(Rejected=sum(Rejected), N=sum(N), Power=Rejected/N, Power.SE=sqrt(Power*(1-Power))/sqrt(N))
        plot.data <- "uniform"
        df <- summary %>%
            filter(Method %in% method.values) %>%
            mutate(Method = factor(Method, method.values, method.labels)) %>%
            mutate(Dimensions = sprintf("m=%d, n=%d", m, n), Dimensions=factor(Dimensions, dimensions.values)) %>%            
            mutate(Alterantive = sprintf("Lehmann's alternative, k=%d", k-1))
        pp <- df %>%
            ggplot(aes(x=m, y=Power, color=Method, shape=Method)) +
            geom_point(alpha=0.75) +
            geom_line(alpha=0.75) +
            ##geom_errorbar(aes(ymin=Power-2*Power.SE, ymax=Power+2*Power.SE), alpha=0.5) +
            geom_hline(yintercept=0.1, linetype=2) +
            facet_grid(.~Alterantive, scales="free") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            ylim(0,NA) +
            #xlim(NA,0.15) +
            xlab("Calibration sample size") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        plot.file.1 <- sprintf("figures/lehmann_%s_cal.pdf", plot.data)
        ggsave(pp, file=plot.file.1, height=3.5, width=7, units="in")
        ## Make table
        for(short.table in c(TRUE,FALSE)) {
            if(short.table) {
                m.values <- c(100,500,1000)
            } else {
                m.values <- seq(100,1000,by=100)
            }
            df.tab <- df %>%
                filter(m %in% m.values) %>%
                mutate(`Calibration size`=m) %>% select(-m) %>%
                select(Dimensions, k, `Calibration size`, Method, Power, Power.SE) %>%
                group_by(Dimensions, k, `Calibration size`) %>%
                mutate(Max=max(Power)) %>%
                ungroup() %>%
                mutate(Power.str = sprintf("%.4f (%.4f)", Power, 2*Power.SE),
                       Power.str = ifelse(Power==Max, sprintf("\\textbf{%s}", Power.str), Power.str)) %>%
                select(k, `Calibration size`, Method, Power.str) %>%       
                pivot_wider(names_from = c(Method), values_from = Power.str, names_sort=TRUE) %>%
                arrange(k, `Calibration size`)           
            tab <- df.tab %>%
                select(-k) %>%
                kable(format="latex", booktabs=TRUE, align = 'c', escape=FALSE) %>%                
                pack_rows("Lehmann's alternative with k=1", start_row = 1, end_row = cumsum(table(df.tab$k))[1]) %>%
                pack_rows("Lehmann's alternative with k=2", start_row = cumsum(table(df.tab$k))[1]+1, end_row = cumsum(table(df.tab$k))[2]) %>%
                pack_rows("Lehmann's alternative with k=3", start_row = cumsum(table(df.tab$k))[2]+1, end_row = cumsum(table(df.tab$k))[3]) %>%
                column_spec(c(1), width = "5em") %>%
                add_header_above(c(" "=1, "Testing procedure" = 6))                
            writeLines(tab, sprintf("tables/lehmann_short%s_cal.tex", short.table))
        }
    }

    results <- load_data_la(92)
    plot.lehmann.new.calibration(reload=FALSE)

    
}

if(plot.lehmann) {

    make_plot_lower_bound_lehmann <- function(plot.data, plot.n_cal, plot.n_test, reload=FALSE, save.plot=TRUE) {
        init_settings(idx.exclude=NULL)
        if(reload){
            results <- load_data(1000)
        }
        summary <- results %>%
            filter(Data==plot.data, n_cal==plot.n_cal, n_test==plot.n_test) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, selected_num, k, n_cal, n_test, prop_out, n_out, Alpha, method) %>%
            summarise(LB=quantile(lower_bound, 0.9), LB.se=sd(lower_bound)/sqrt(n()), N=n())
        df <- summary %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels))
        pp <- df %>%
            ggplot(aes(x=n_out, y=LB, color=Method, shape=Method)) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se)), width=0.1, alpha=0.5) +
            geom_abline(slope=1, intercept=0, linetype=2) +
            facet_grid(k~selected_num) +#, labeller="label_botval") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("True number of outliers") +
            ylab("90% lower bound") +
            scale_x_continuous(trans='log10') +
            scale_y_continuous(trans='log10') +
                                        #            xlim(NA,0.5) +
                                        #            ylim(NA,0.5) +
            theme_bw()
        if(save.plot) {
            plot.file.1 <- sprintf("figures/lehmann_lb_%s_n%d_n%d.png", plot.data, plot.n_cal, plot.n_test)
            ggsave(pp, file=plot.file.1, height=4, width=7, units="in")
        } else {
            print(pp)
        }
    }

    make_plot_power_lehmann <- function(plot.data, plot.n_cal, plot.n_test, reload=FALSE, save.plot=TRUE) {
        init_settings(idx.exclude=NULL)
        if(reload){
            results <- load_data(1000)
        }
        summary <- results %>%
            filter(Data==plot.data, n_cal==plot.n_cal, n_test==plot.n_test) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, selected_num, k, n_cal, n_test, prop_out, n_out, Alpha, method) %>%
            summarise(Power=mean(lower_bound>0), SE=sd(lower_bound>0)/sqrt(n()), N=n())
        df <- summary %>%
            filter(method %in% method.values) %>%
            mutate(Method = factor(method, method.values, method.labels))
        pp <- df %>%
            ggplot(aes(x=n_out, y=Power, color=Method, shape=Method)) +
            geom_point() +
            geom_line() +
#            geom_errorbar(aes(ymin=(Power-2*SE), ymax=(Power+2*SE)), width=0.1, alpha=0.5) +
            geom_abline(slope=0, intercept=0.1, linetype=2) +
            facet_grid(k~selected_num) +#, labeller="label_botval") +
            scale_color_manual(values=color.scale) +
            scale_shape_manual(values=shape.scale) +
            scale_alpha_manual(values=alpha.scale) +
            xlab("True number of outliers") +
            ylab("90% lower bound") +
#            scale_x_continuous(trans='log10') +
#            scale_y_continuous(trans='log10') +
                                        #            xlim(NA,0.5) +
                                        #            ylim(NA,0.5) +
            theme_bw()
        if(save.plot) {
            plot.file.1 <- sprintf("figures/lehmann_power_%s_n%d_n%d.png", plot.data, plot.n_cal, plot.n_test)
            ggsave(pp, file=plot.file.1, height=4, width=7, units="in")
        } else {
            print(pp)
        }
    }
    
    results <- load_data(1000)

    for(plot.data in c("uniform", "exponential", "normal")) {
        for(n.plot in c(1000)) {
            ##make_plot_lower_bound_lehmann(plot.data, plot.n_cal=n.plot, plot.n_test=n.plot, reload=FALSE, save.plot=TRUE)
            make_plot_power_lehmann(plot.data, plot.n_cal=n.plot, plot.n_test=n.plot, reload=FALSE, save.plot=TRUE)
        }
    }
    
}


if(plot.lhco_selection) {

    make_plot_lower_bound_lhco_sel <- function(setup, plot.n_train, reload=FALSE, plot.quantile=0.5, plot.simple=FALSE, plot.naive=FALSE) {
        if(plot.simple) {
            init_settings(idx.exclude=c(2,4,6))
        } else {
            init_settings(idx.exclude=c(6))
        }
        if(reload) {
            results <- load_data(setup)
        }
        summary <- results %>%
            filter(n_train == plot.n_train) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3"), names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier, method) %>%
            summarise(LB=quantile(lower_bound, plot.quantile), LB.se=sd(lower_bound)/sqrt(n()), N=n())
        if(plot.naive) {
            summary.naive <- results %>%
                filter(n_train == plot.n_train) %>%
                pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3"), names_to="method", values_to="lower_bound") %>%
                group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, seed, Seed, Repetition) %>%
                summarise(Method="greedy", method="greedy", Classifier="greedy", lower_bound=max(lower_bound)) %>%
                group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier, method) %>%
                summarise(LB=quantile(lower_bound, plot.quantile), LB.se=sd(lower_bound)/sqrt(n()), N=n())
            summary <- rbind(summary, summary.naive)
            method.values.tmp <- c(method.values, "greedy")
            method.labels.tmp <- c(method.labels, "Naive heuristic")
            color.scale.tmp <- c(color.scale, cbPalette[2])
            shape.scale.tmp <- c(shape.scale, 4)
            alpha.scale.tmp <- c(alpha.scale, 0.5)
        } else {
            method.values.tmp <- method.values
            method.labels.tmp <- method.labels
            color.scale.tmp <- color.scale
            shape.scale.tmp <- shape.scale
            alpha.scale.tmp <- alpha.scale
        }
        prop.out.values <- c(0.05,0.1,0.15)
        N_out.values <- paste(10000*prop.out.values, "outliers")
        df <- summary %>%
            filter(prop_out %in% prop.out.values) %>%
            filter(method %in% method.values.tmp) %>%
            mutate(Method = factor(method, method.values.tmp, method.labels.tmp)) %>%
            mutate(Classifier = factor(Classifier, classifier.values, classifier.labels)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
            mutate(N_out = factor(N_out, N_out.values, N_out.values)) %>%
            mutate(Data = factor(Data, data.values, data.labels))
        df.ref <- results %>%
            filter(n_train == plot.n_train) %>%
            filter(prop_out %in% prop.out.values) %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier) %>%
            summarise(LB=quantile(n_out_sel, plot.quantile)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
            mutate(N_out = factor(N_out, N_out.values, N_out.values)) %>%
            mutate(Data = factor(Data, data.values, data.labels))
        df.ghost <- tibble(prop_out = c(0.05,0.05), selected_num=c(0,0), n_test=c(10000,10000), LB=c(0,1)) %>%
            mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
            mutate(N_out = factor(N_out, N_out.values, N_out.values))
        x.min <- min(df$selected_num/df$n_test)
        x.max <- max(df$selected_num/df$n_test)
        pp <- df %>%
            ggplot(aes(x=selected_num/n_test, y=LB)) +
            geom_point(aes(color=Method, shape=Method, alpha=Method)) +
            geom_line(aes(color=Method, alpha=Method)) +
            geom_line(data=df.ref, linetype=2, aes(x=selected_num/n_test, y=LB), color="black", alpha=1) +
            geom_point(data=df.ghost, aes(x=selected_num/n_test, y=LB), alpha=0) +
            ##geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se), color=Method), width=0.01, alpha=0.5) +
            facet_wrap(.~N_out, labeller="label_value", scale="free", nrow=1) +
            scale_color_manual(values=color.scale.tmp) +
            scale_shape_manual(values=shape.scale.tmp) +
            scale_alpha_manual(values=alpha.scale.tmp) +
            xlab("Proportion of top scores in selected test set") +
            ylab("90% lower bound") +
            scale_x_continuous(trans='log10', limits=c(0.01,1)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
            guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))          
        plot.file.1 <- sprintf("figures/lhco%d_nt%d_lower_bound_sel_simple%s_q%s.pdf", setup, plot.n_train, plot.simple, plot.quantile)
        ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
    }


    results <- load_data(1100)

    make_plot_lower_bound_lhco_sel(1100, plot.n_train=10000, plot.simple=FALSE, plot.quantile=0.5, reload=FALSE)
    make_plot_lower_bound_lhco_sel(1100, plot.n_train=10000, plot.simple=FALSE, plot.quantile=0.9, reload=FALSE)
    make_plot_lower_bound_lhco_sel(1100, plot.n_train=100000, plot.simple=FALSE, plot.quantile=0.5, reload=FALSE)
    make_plot_lower_bound_lhco_sel(1100, plot.n_train=100000, plot.simple=FALSE, plot.quantile=0.9, reload=FALSE)

    ##make_plot_lower_bound_lhco_sel(1100, plot.n_train=100000, reload=FALSE)

}


if(plot.lhco_2) {

    make_plot_lower_bound_lhco_2 <- function(reload=FALSE, plot.n_train, plot.n_cal, include.BH=FALSE,  plot.classifier="auto",
                                             plot.simple=FALSE, plot.naive=FALSE, save.plot=True) {
        init_settings_local <- function(idx.exclude=NULL) {
            cbPalette <<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#6e57d2", "red")
            method.values <<- c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4", "lb_auto")
            method.labels <<- c("Simes", "Storey-Simes", "Fisher", "Mann-Whitney-Wilcoxon", "WMW (k=2)", "WMW (k=3)", "Adaptively selected")
            classifier.values <<- c("occ-auto", "bc-auto", "auto")
            classifier.labels <<- c("One-Class", "Binary", "Automatic")
            color.scale <<- cbPalette[c(1,1,4,3,6,7,8)]
            shape.scale <<- c(2,6,3,1,0,9,8)
            alpha.scale <<- c(0.75,0.75,0.75,0.75,0.75,0.75,1)
            data.values <<- c("pendigits", "creditcard", "cover", "shuttle", "mammography", "aloi")
            data.labels <<- c("Pendigits", "Creditcard", "Covertype", "Shuttle", "Mammography", "ALOI")
            if(length(idx.exclude)>0) {
                ## Exclude Storey-simes
                method.values <<- method.values[-idx.exclude]
                method.labels <<- method.labels[-idx.exclude]
                color.scale <<- color.scale[-idx.exclude]
                shape.scale <<- shape.scale[-idx.exclude]
                alpha.scale <<- alpha.scale[-idx.exclude]
            }
        }
        if(plot.simple) {
            init_settings_local(idx.exclude=c(2,4,5,6))
        } else {
            init_settings_local(idx.exclude=NULL)            
        }
        if(reload){
            results <- load_data(100)
        }
        pow.str <- "Power (global null, no outliers)"
        ##        key.values <- c("Median", "Quantile.90", pow.str)
        key.values <- c("Median", pow.str)
        ##        key.labels <- c("Lower bound (median)", "Lower bound (90th quant.)", pow.str)
        key.labels <- c("Lower bound for num. outliers", pow.str)
        disc.str <- "Discoveries (10% FDR)"
        df.lb <- results %>%
            filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier==plot.classifier) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method, tune_size, selection) %>%
            summarise(Median=median(lower_bound), Quantile.90=quantile(lower_bound, 0.9), SE=sd(lower_bound)/sqrt(n()), N=n()) %>%
            pivot_longer(c("Median", "Quantile.90"), names_to="Key", values_to="Value")
        if(plot.naive) {
            method.values.tmp <- c(method.values, "greedy")
            method.labels.tmp <- c(method.labels, "Cherry picking")
            color.scale.tmp <- c(color.scale, cbPalette[2])
            shape.scale.tmp <- c(shape.scale, 4)
            alpha.scale.tmp <- c(alpha.scale, 0.5)
            df.lb.greedy <- results %>%
                filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier!="auto") %>%
                pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
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
            color.scale.tmp <- color.scale
            shape.scale.tmp <- shape.scale
            alpha.scale.tmp <- alpha.scale
            df.lb <- df.lb %>%
                filter(method %in% method.values) %>%
                mutate(Method = factor(method, method.values, method.labels), Target="Closed testing")
        }
        df.pow <- results %>%
            filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier==plot.classifier) %>%
            pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
                         names_to="method", values_to="lower_bound") %>%
            group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method, tune_size, selection) %>%
            summarise(Key=pow.str, Value=mean(lower_bound>0), SE=sd(lower_bound>0)/sqrt(n()), N=n())
        if(plot.naive) {
            df.pow.greedy <- results %>%
                filter(n_train == plot.n_train, n_cal == plot.n_cal, Classifier!="auto") %>%
                pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw_k2", "lb_wmw_k3", "lb_wmw_k4"),
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
        df.range <- tibble(Key=c(pow.str, pow.str, "Median", "Median", "Quantile.90", "Quantile.90"), Value=c(0,1,0,1500,0,1500),
                           n_out=c(0,1500,0,1500,0,1500), method="lb_auto") %>%
            mutate(Method = factor(method, method.values.tmp, method.labels.tmp), Target="Closed testing") %>%
            filter(Key %in% key.values) %>%
            mutate(Key = factor(Key, key.values, key.labels))
        df.ref <- tibble(Key=c(pow.str, pow.str, "Median", "Median", "Quantile.90", "Quantile.90"),
                         Value=c(0.1,0.1,0,1500,0,1500), n_out=c(0,1500,0,1500,0,1500), method="lb_auto", Target="Closed testing") %>%
            mutate(Method = factor(method, method.values.tmp, method.labels.tmp)) %>%
            filter(Key %in% key.values) %>%
            mutate(Key = factor(Key, key.values, key.labels))
        if(include.BH) {
            pp <- df %>%
                ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method, linetype=Target)) +
                geom_point() +
                geom_line() +
                ##geom_errorbar(aes(ymin=(Value-2*SE), ymax=(Value+2*SE)), width=0.1, alpha=0.5) +
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
                guides(linetype = "none", color=guide_legend(title="Aggregation method"), shape=guide_legend(title="Aggregation method"), alpha=guide_legend(title="Aggregation method"))
        } else {
            pp <- df %>%
                ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method, linetype=Target)) +
                geom_point() +
                geom_line() +
                ##geom_errorbar(aes(ymin=(Value-2*SE), ymax=(Value+2*SE)), width=0.1, alpha=0.5) +
                geom_point(data=df.range, aes(x=n_out, y=Value, color=Method, shape=Method), alpha=0) +
                geom_line(data=df.ref, aes(x=n_out, y=Value), linetype=2, color="black", alpha=1) +
                                        #                geom_line(data=df.fdr, aes(x=n_out, y=Value, linetype=Target), color="black", alpha=1) +
                facet_wrap(.~Key, labeller="label_value", scales="free") +
                scale_color_manual(values=color.scale.tmp) +
                scale_shape_manual(values=shape.scale.tmp) +
                scale_alpha_manual(values=alpha.scale.tmp) +
                scale_linetype_manual(values=c(1,3)) +
                xlab("True number of outliers") +
                ylab("") +
                theme_bw() +
                guides(linetype = "none", color=guide_legend(title="Conformal testing approach"),
                       shape=guide_legend(title="Conformal testing approach"),
                       alpha=guide_legend(title="Conformal testing approach"))
        }
        if(save.plot) {
            plot.file.1 <- sprintf("figures/demo_outliers.pdf", plot.n_train, include.BH, plot.simple)
            ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
        } else {
            print(pp)
        }
    }

    results <- load_data(100)


    make_plot_lower_bound_lhco_2(reload=FALSE, plot.n_train=10000, plot.n_cal=2000, include.BH=TRUE, plot.classifier="auto",
                                 plot.simple=TRUE, plot.naive=FALSE, save.plot=TRUE)


    #make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=100000, plot.n_cal=2000, include.BH=FALSE, save.plot=TRUE)
    #make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=10000, plot.n_cal=2000, include.BH=TRUE, save.plot=TRUE)
    #make_plot_lower_bound_lhco(reload=FALSE, plot.n_train=100000, plot.n_cal=2000, include.BH=TRUE, save.plot=TRUE)



}
