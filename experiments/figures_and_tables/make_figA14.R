## -----------------------------------------------------------------------------
## Required results : figA14
## Produces         : Figure A14
## -----------------------------------------------------------------------------

source("plotting_utils_A.R")

## ---- Plot function: lower bound vs proportion of top scores selected -------
make_plot_lower_bound_lhco_sel <- function(fig.name, plot.n_train, plot.quantile=0.5,
                                           plot.simple=FALSE, plot.naive=FALSE, save.plot=TRUE) {
    present_methods = c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw")
    init_settings(present_methods = present_methods)
    summary <- results %>%
        filter(n_train == plot.n_train) %>%
        pivot_longer(all_of(present_methods), names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier, method) %>%
        summarise(LB=quantile(lower_bound, plot.quantile), LB.se=sd(lower_bound)/sqrt(n()), N=n())
    if(plot.naive) {
        summary.naive <- results %>%
            filter(n_train == plot.n_train) %>%
            pivot_longer(all_of(present_methods), names_to="method", values_to="lower_bound") %>%
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
    prop.out.values <- c(0.1, 0.15, 0.25)
    N_out.values <- paste(2000*prop.out.values, "outliers")
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
        summarise(LB=mean(n_out_sel)) %>%
        mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
        mutate(N_out = factor(N_out, N_out.values, N_out.values)) %>%
        mutate(Data = factor(Data, data.values, data.labels))
    df.ghost <- tibble(prop_out = c(0.1, 0.1), selected_num=c(0,0), n_test=c(2000,2000), LB=c(0,1)) %>%
        mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
        mutate(N_out = factor(N_out, N_out.values, N_out.values))
    pp <- df %>%
        ggplot(aes(x=selected_num/n_test, y=LB)) +
        geom_point(aes(color=Method, shape=Method, alpha=Method)) +
        geom_line(aes(color=Method, alpha=Method)) +
        geom_line(data=df.ref, linetype=2, aes(x=selected_num/n_test, y=LB), color="black", alpha=1) +
        geom_point(data=df.ghost, aes(x=selected_num/n_test, y=LB), alpha=0) +
        facet_wrap(.~N_out, labeller="label_value", scale="free", nrow=1) +
        scale_color_manual(values=color.scale.tmp) +
        scale_shape_manual(values=shape.scale.tmp) +
        scale_alpha_manual(values=alpha.scale.tmp) +
        xlab("Proportion of top scores in selected test set") +
        ylab("90% lower bound") +
        scale_x_continuous(trans='log10', limits=c(0.01,1)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        guides(linetype = "none",
               color=guide_legend(title="Local tests"),
               shape=guide_legend(title="Local tests"),
               alpha=guide_legend(title="Local tests"))
    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_lhco_nt%d_lower_bound_sel_simple%s_q%s.pdf",
                               fig.name, plot.n_train, plot.simple, plot.quantile)
        ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
    } else {
        print(pp)
    }
}

results <- load_data("figA14")

## ---- Run --------------------------------------------------------------------
make_plot_lower_bound_lhco_sel(fig.name="figA14", plot.n_train=10000,  plot.quantile=0.5)
