## -----------------------------------------------------------------------------
## Required results : figA11
## Produces         : Figure A11
## -----------------------------------------------------------------------------

source("utils_plotting_B.R")
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---- Plot function: power vs number of outliers, faceted by alternative ----
plot_power <- function(fig.name, save.plot=TRUE) {
    init_settings_power()

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

    ## Significance level
    alpha <- 0.05

    ## Calculate power for different methods and prop_out values
    power_results.raw <- results %>%
        group_by(n_cal, n_test, alternative, Method, prop_out) %>%
        summarize(Power = mean(p.value < alpha),
                  SE = sqrt((Power * (1 - Power)) / n()))

    power_results <- power_results.raw %>%
        filter(Method %in% method.values, alternative %in% alternative.values) %>%
        mutate(Alternative = factor(alternative, alternative.values, alternative.labels),
               Method = factor(Method, method.values, method.labels))

    df <- power_results %>%
        filter(n_cal==n_cal.plot, n_test==n_test.plot, prop_out<=0.25)

    pp <- df %>%
        mutate(n_out = round(prop_out*n_test)) %>%
        ggplot(aes(x = n_out, y = Power, color = Method, shape = Method)) +
        geom_line() +
        geom_point() +
        geom_hline(yintercept = alpha, linetype = 2) +
        facet_wrap(.~Alternative, nrow=2, labeller="label_value") +
        ylim(0,1) +
        theme_bw(base_size = 15) +
        scale_color_manual(values = color.by.label) +
        scale_shape_manual(values = shape.by.label) +
        labs(x = "Number of Outliers",
             y = "Power",
             color = "Method") +
        theme(legend.position = "bottom")

    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_power_ncal%d_ntest%d.pdf",
                               fig.name, n_cal.plot, n_test.plot)
        ggsave(filename = plot.file.1, plot = pp, width = 10, height = 5)
    } else {
        print(pp)
    }
}

results <- load_data("figA11")

## ---- Run --------------------------------------------------------------------
plot_power(fig.name="figA11")
