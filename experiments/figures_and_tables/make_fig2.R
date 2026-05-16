## -----------------------------------------------------------------------------
## Required results : fig2
## Produces         : Figure 2, Figure A4, Figure A7
##
## Required results : figA5
## Produces         : Figure A5, Figure A9
## -----------------------------------------------------------------------------

source("plotting_utils_A.R")

## ---- Plot function: 90% lower bound vs. true number of outliers, faceted by classifier ----
make_plot_lower_bound_proportion <- function(fig.name, plot.quantile=0.5, include_BH=FALSE) {
    if (include_BH) {
        present_methods <- c("disc_bh", "lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_wmw")
    } else {
        present_methods <- c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_higher_k2")
    }
    init_settings(present_methods = present_methods)
    summary <- results %>%
        pivot_longer(all_of(present_methods),
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
    if (include_BH) {
        plot.file.1 <- sprintf("figures/%s_lower_bound_q%s_bh.pdf", fig.name, plot.quantile)
    } else {
        plot.file.1 <- sprintf("figures/%s_lower_bound_q%s.pdf", fig.name, plot.quantile)
    }
    ggsave(pp, file=plot.file.1, height=2.25, width=7.5, units="in")
}

## ---- Run --------------------------------------------------------------------

## Setup 1: circles-mixed synthetic data
results <- load_data("fig2")

## Paper Figure 2: median lower bound, no BH overlay
make_plot_lower_bound_proportion(fig.name="fig2", plot.quantile=0.5)

## Appendix: with BH overlay
make_plot_lower_bound_proportion(fig.name="figA4", plot.quantile=0.5, include_BH=TRUE)

## Appendix: 90th-quantile version
make_plot_lower_bound_proportion(fig.name="figA7", plot.quantile=0.9)

## Setup 2: binomial synthetic data (appendix)
results <- load_data("figA5")

make_plot_lower_bound_proportion(fig.name="figA5", plot.quantile=0.5)
make_plot_lower_bound_proportion(fig.name="figA9", plot.quantile=0.9)
