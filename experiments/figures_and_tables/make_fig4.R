## -----------------------------------------------------------------------------
## Required results : fig4
## Produces         : Figure 4
## -----------------------------------------------------------------------------

source("plotting_utils_A.R")

## ---- Plot function: 90% lower bound vs. true number of outliers, faceted by classifier ----
## Used for paper Figure 4 (adversarial-anomaly experiments).
make_plot_lower_bound_adversary <- function(fig.name, plot.quantile=0.5) {
    init_settings(idx.exclude=NULL)
    summary <- results %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_higher_k2", "lb_lmp"),
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
        ##geom_errorbar(aes(ymin=(LB-2*LB.se), ymax=(LB+2*LB.se)), width=0.01, alpha=0.5) +
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
    plot.file.1 <- sprintf("figures/%s_lower_bound_q%s.pdf", fig.name, plot.quantile)
    ggsave(pp, file=plot.file.1, height=2.25, width=7, units="in")
}

## ---- Plot function: combined Lower Bound + Power panels, faceted by classifier ----
## Companion figure with both metrics in one chart.
make_combined_plot_adversary <- function(fig.name, plot.quantile=0.5, plot.alpha=0.1) {
    init_settings(idx.exclude=NULL)
    ## Lower-bound summary
    lb_summary <- results %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_higher_k2", "lb_lmp"),
                     names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
        summarise(Value=quantile(lower_bound, plot.quantile), Value.se=sd(lower_bound)/sqrt(n()), .groups = 'drop') %>%
        mutate(metric = "Lower Bound")
    ## Power summary (uses pval_* columns, renamed to lb_* to match method.values)
    power_summary <- results %>%
        pivot_longer(c("pval_simes", "pval_storey_simes", "pval_fisher", "pval_auto", "pval_wmw", "pval_higher_k2", "pval_lmp"),
                     names_to="method", values_to="pval") %>%
        mutate(method = str_replace(method, "^pval_", "lb_")) %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
        summarise(Value=mean(pval<=plot.alpha), Value.se=sd(pval<=plot.alpha)/sqrt(n()), .groups = 'drop') %>%
        mutate(metric = "Power (global null)")
    ## Combine
    combined_summary <- bind_rows(lb_summary, power_summary)
    df <- combined_summary %>%
        filter(method %in% method.values) %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Classifier = factor(Classifier, classifier.values, classifier.labels))
    x.max <- max(df$n_out)
    ## Reference lines and axis-pinning ghost points, per panel
    data.ref.pow <- tibble(n_out=c(0, max(df$n_out)), metric="Power (global null)", Value=c(0.1,0.1), Method="Adaptive") %>%
        mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound")))
    data.ref.lb <- tibble(n_out=c(0, max(df$n_out)), metric="Lower Bound", Value=c(0, max(df$n_out)), Method="Adaptive") %>%
        mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound")))
    data.lim.pow <- tibble(n_out=c(0, 0), metric="Power (global null)", Value=c(0,1), Method="Adaptive") %>%
        mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound")))
    data.lim.lb <- tibble(n_out=c(0, 0), metric="Lower Bound", Value=c(0, max(df$n_out)), Method="Adaptive") %>%
        mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound")))
    pp <- df %>%
        mutate(metric = factor(metric, levels=c("Power (global null)", "Lower Bound"))) %>%
        ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method)) +
        geom_point() +
        geom_line() +
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
        theme_bw() +
        guides(linetype = "none", color=guide_legend(title="Local tests"), shape=guide_legend(title="Local tests"), alpha=guide_legend(title="Local tests"))
    plot.file.1 <- sprintf("figures/%s_combined_plot_q%s_alpha%s.pdf", fig.name, plot.quantile, plot.alpha)
    ggsave(pp, file=plot.file.1, height=3.5, width=7, units="in")
}

## ---- Run --------------------------------------------------------------------

results <- load_data("fig4")

## Paper Figure 4: 90th-percentile lower bound under adversarial anomalies
make_plot_lower_bound_adversary(fig.name="fig4", plot.quantile=0.9)

## Companion: combined Lower Bound + Power panels
##make_combined_plot_adversary(fig.name="fig4", plot.quantile=0.9)
