## -----------------------------------------------------------------------------
## Required results : figA12
## Produces         : Figure A12
## -----------------------------------------------------------------------------

source("utils_plotting_B.R")
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---- Plot function: lower bound vs number of outliers, faceted by alternative
plot_lb <- function(fig.name, save.plot=TRUE) {
    init_settings_lb()

    n_cal.plot <- 500
    n_test.plot <- 200
    alternative.values <- c("uniform", "lehmann_k2",
                            "beta_0.25_0.25", "beta_10_10",
                            "normal_1.5_1", "normal_-1.5_1",
                            "normal_0_2", "normal_0_0.25")
    alternative.labels <- c("Uniform (null)", "Lehmann (k=3)",
                            "Beta (overdispersed)", "Beta (underdispersed)",
                            "Normal (positive shift)", "Normal (negative shift)",
                            "Normal (overdispersed)", "Normal (underdispersed)")

    ## Aggregate lower bound across replicates
    lb_results.raw <- results %>%
        mutate(n.out = round(n_test*prop_out)) %>%
        group_by(n_cal, n_test, alternative, Method, prop_out, n.out) %>%
        summarize(LB = median(Lower),
                  SE = sd(Lower)/sqrt(n())) %>%
        mutate(SE = ifelse(is.na(SE), 0, SE))

    lb_results <- lb_results.raw %>%
        filter(n_test == n_test.plot) %>%
        filter(Method %in% method.values, alternative %in% alternative.values) %>%
        mutate(Alternative = factor(alternative, alternative.values, alternative.labels),
               Method = factor(Method, method.values, method.labels))

    df <- lb_results %>%
        filter(n_cal==n_cal.plot, n_test==n_test.plot)

    pp <- df %>%
        ggplot(aes(x = n.out, y = LB, color = Method, shape = Method)) +
        geom_line() +
        geom_point() +
        geom_abline(slope=1, linetype = 2) +
        facet_wrap(.~Alternative, nrow=2, labeller="label_value") +
        theme_bw(base_size = 15) +
        scale_color_manual(values = color.by.label) +
        scale_shape_manual(values = shape.by.label) +
        labs(x = "Number of Outliers",
             y = "Lower bound",
             color = "Method") +
        theme(legend.position = "bottom")

    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_lb_ncal%d_ntest%d.pdf",
                               fig.name, n_cal.plot, n_test.plot)
        ggsave(filename = plot.file.1, plot = pp, width = 10, height = 5)
    } else {
        print(pp)
    }
}

results <- load_data("figA12")

## ---- Run --------------------------------------------------------------------
plot_lb(fig.name="figA12")
