## -----------------------------------------------------------------------------
## Required results : figA10
## Produces         : Figure A10
## -----------------------------------------------------------------------------

source("plotting_utils_A.R")

## ---- Plot function: lower bound vs mixture proportion, faceted by classifier
make_plot_lower_bound_mixture <- function(fig.name, plot.p, plot.n_train, plot.signal, save.plot=TRUE) {
    init_settings(idx.exclude=c(6))
    summary <- results %>%
        separate(Data, into = c("Data","Mixture"), sep = "-") %>%
        mutate(Mixture=parse_number(Mixture)) %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_higher_k2"),
                     names_to="method", values_to="lower_bound") %>%
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
        facet_grid(.~Classifier, labeller="label_both") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_alpha_manual(values=alpha.scale) +
        xlab("Proportion of data points from binomial distribution") +
        ylab("90% lower bound") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        guides(linetype = "none",
               color=guide_legend(title="Local tests"),
               shape=guide_legend(title="Local tests"),
               alpha=guide_legend(title="Local tests"))
    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_synthetic3_p%d_n%d_s%.2f_lower_bound.pdf",
                               fig.name, plot.p, plot.n_train, plot.signal)
        ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
    } else {
        print(pp)
    }
}

results <- load_data("figA10")

## ---- Run --------------------------------------------------------------------
make_plot_lower_bound_mixture(fig.name="figA10", plot.p=100, plot.n_train=1000, plot.signal=3)
