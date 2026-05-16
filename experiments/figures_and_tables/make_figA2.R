## -----------------------------------------------------------------------------
## Required results : figA2
## Produces         : Figure A2
## -----------------------------------------------------------------------------

source("plotting_utils_A.R")

## ---- Plot function: tuning curves faceted by classifier x #outliers --------
make_plot_lower_bound_tuning <- function(fig.name, plot.quantile=0.5, save.plot=TRUE) {
    init_settings(idx.exclude=c(6))
    summary <- results %>%
        mutate(tune_size,
               `True number of outliers`=n_out) %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_higher_k2"),
                     names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, tune_size, prop_out,
                 `True number of outliers`, Alpha, Classifier, method) %>%
        summarise(LB=quantile(lower_bound, plot.quantile),
                  LB.se=sd(lower_bound)/sqrt(n()))
    df <- summary %>%
        filter(method %in% method.values) %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
    pp <- df %>%
        ggplot(aes(x=tune_size, y=LB, color=Method, shape=Method, alpha=Method)) +
        geom_point() +
        geom_line() +
        geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se)),
                      width=0.01, alpha=0.5) +
        facet_grid(Classifier~`True number of outliers`, labeller="label_both") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_alpha_manual(values=alpha.scale) +
        xlab("Proportion of calibration samples used for tuning") +
        ylab("90% lower bound") +
        scale_x_continuous(trans='log10', limits=c(0.01,1)) +
        theme_bw() +
        guides(linetype = "none",
               color=guide_legend(title="Local tests"),
               shape=guide_legend(title="Local tests"),
               alpha=guide_legend(title="Local tests")) +
        theme(panel.spacing.x = unit(3, "mm"))
    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_synthetic1t_lower_bound_q%s_tuning.pdf",
                               fig.name, plot.quantile)
        ggsave(pp, file=plot.file.1, height=5, width=8.5, units="in")
    } else {
        print(pp)
    }
}

results <- load_data("figA2")

## ---- Run --------------------------------------------------------------------
make_plot_lower_bound_tuning(fig.name="figA2", plot.quantile=0.5)
