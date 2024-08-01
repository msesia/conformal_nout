options(width = 300)

library(tidyverse)
#library(latex2exp)
library(RColorBrewer)
#library(kableExtra)

plot.power_1 <- TRUE
plot.lb_1 <- TRUE

load_data <- function(setup) {
    idir <- sprintf("results_hpc/setup%d", setup)
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}

init_settings <- function(idx.exclude=NULL) { 
    method.values <- c("Fisher", "WMW", "Shirashi_oracle", "Shirashi_ghat_betamix", "Shirashi_ghat_betamix_inc")
    method.labels <- c("Fisher", "WMW", "LMP (oracle)", "LMP (empirical)", "LMP (empirical, monotone)")
    alternative.values <- c("uniform", "lehmann_k2", "beta_0.5_0.5", "beta_4_4", "normal_0.5_1", "normal_-0.5_1", "normal_0_1.5", "normal_0_0.5")
    alternative.labels <- c("Uniform (null)", "Lehmann", "Beta (overdispersed)", "Beta (underdispersed)", "Normal (positive shift)", "Normal (negative shift)",
                            "Normal (overdispersed)", "Normal (underdispersed)")
    ## Manual color and shape scales
    colors <- c("Fisher" = "#66A9D2",
                "WMW" = "#E69F00", 
                "LMP (oracle)" = "#00441B",    # Very dark green
                "LMP (empirical)" = "#238B45", # Darker green
                "LMP (empirical, monotone)" = "#41AB5D") # Dark green
    shapes <- c("Fisher" = 1, "WMW" = 2, "LMP (oracle)" = 8, 
                "LMP (empirical)" = 15, "LMP (empirical, monotone)" = 16)
}

plot_lb_1 <- function() {

    results <- load_data(2)

    method.values <- c("Fisher", "WMW", "Shirashi_oracle", "Shirashi_ghat_betamix")
    method.labels <- c("Fisher", "WMW", "LMP (oracle)", "LMP (empirical, monotone)")
##    alternative.values <- c("uniform", "lehmann_k2", "beta_0.5_0.5", "beta_4_4", "normal_0.5_1", "normal_-0.5_1", "normal_0_0.5", "normal_0_1.5")
    alternative.values <- c("uniform", "lehmann_k2", "beta_0.25_0.25", "beta_4_4", "normal_2_1", "normal_-2_1", "normal_0_2", "normal_0_0.25")
    alternative.labels <- c("Uniform (null)", "Lehmann", "Beta (overdispersed)", "Beta (underdispersed)", "Normal (positive shift)", "Normal (negative shift)",
                            "Normal (overdispersed)", "Normal (underdispersed)")
    ## Manual color and shape scales
    colors <- c("Fisher" = "#66A9D2",
                "WMW" = "#E69F00", 
                "LMP (oracle)" = "#00441B",    # Very dark green
                "LMP (empirical, monotone)" = "#238B45") # Darker green
    shapes <- c("Fisher" = 1, "WMW" = 2, "LMP (oracle)" = 8, 
                "LMP (empirical)" = 15, "LMP (empirical, monotone)" = 16)
    
    ## Calculate power for different methods and prop_out values
    lb_results.raw <- results %>%
        mutate(n.out = round(n_test*prop_out)) |>
        group_by(n_cal, n_test, alternative, Method, prop_out, n.out) %>%
        summarize(
            LB = mean(Lower),
            SE = sd(Lower)/sqrt(n())
        ) |>
    mutate(SE = ifelse(is.na(SE), 0, SE))

    lb_results <- lb_results.raw |>
    filter(n_test == 200) |>
    filter(Method %in% method.values, alternative %in% alternative.values) |>
    mutate(Alternative = factor(alternative, alternative.values, alternative.labels),
           Method = factor(Method, method.values, method.labels))

    ## Function to plot and save the power plot for a given alternative value
    plot_lb_for_n <- function(n_cal.plot, n_test.plot) {
        ## Filter for the specified alternative
        df <- lb_results |>
        filter(n_cal==n_cal.plot, n_test == n_test.plot)
        ## Make plot
        pp <- df |>
        ggplot(aes(x = n.out, y = LB, color = Method, shape = Method)) +
            geom_line() +
            geom_point() +
            geom_errorbar(aes(ymin = LB - 2*SE, ymax = LB + 2*SE), width = 0.02) +
            geom_abline(slope=1, linetype = 2) +
            facet_wrap(.~Alternative, nrow=2, labeller="label_value") +
            theme_bw(base_size = 15) +
            scale_color_manual(values = colors) +
            scale_shape_manual(values = shapes) +
            labs(#title = paste("LB of Different Methods for Global Testing"),
                 #subtitle = sprintf("Calibration size: %d, Test size: %d", n_cal.plot, n_test.plot),
                 x = "Number of Outliers",
                 y = "Lower bound",
                 color = "Method") +
            theme(legend.position = "bottom")       
        ## Save the plot as a PNG file
        filename <- sprintf("figures/lb_ncal%d_ntest%d.png", n_cal.plot, n_test.plot)
        ggsave(filename = filename, plot = pp, width = 10, height = 5)
    }

    plot_lb_for_n(500,200)
    
}

if(plot.lb_1) {
    plot_lb_1()
}



plot_power_1 <- function() {

    results <- load_data(1)
    

    method.values <- c("Fisher", "WMW", "Shirashi_oracle", "Shirashi_ghat_betamix", "Shirashi_ghat_betamix_inc")
    method.labels <- c("Fisher", "WMW", "LMP (oracle)", "LMP (empirical)", "LMP (empirical, monotone)")
    alternative.values <- c("uniform", "lehmann_k2", "beta_0.25_0.25", "beta_4_4", "normal_2_1", "normal_-2_1", "normal_0_2", "normal_0_0.25")
    alternative.labels <- c("Uniform (null)", "Lehmann", "Beta (overdispersed)", "Beta (underdispersed)", "Normal (positive shift)", "Normal (negative shift)",
                            "Normal (overdispersed)", "Normal (underdispersed)")
    ## Manual color and shape scales
    colors <- c("Fisher" = "#66A9D2",
                "WMW" = "#E69F00", 
                "LMP (oracle)" = "#00441B",    # Very dark green
                "LMP (empirical)" = "#238B45", # Darker green
                "LMP (empirical, monotone)" = "#41AB5D") # Dark green
    shapes <- c("Fisher" = 1, "WMW" = 2, "LMP (oracle)" = 8, 
                "LMP (empirical)" = 15, "LMP (empirical, monotone)" = 16)

    ## Significance level
    alpha <- 0.1

    ## Calculate power for different methods and prop_out values
    power_results.raw <- results %>%
        group_by(n_cal, n_test, alternative, Method, prop_out) %>%
        summarize(
            Power = mean(p.value < alpha),
            SE = sqrt((Power * (1 - Power)) / n())
        )

    
    power_results <- power_results.raw |>
    filter(n_test == 200) |>
    filter(Method %in% method.values, alternative %in% alternative.values) |>
    mutate(Alternative = factor(alternative, alternative.values, alternative.labels),
           Method = factor(Method, method.values, method.labels))

    ## Function to plot and save the power plot for a given alternative value
    plot_power_for_n <- function(n_cal.plot, n_test.plot) {
        ## Filter for the specified alternative
        df <- power_results |>
        filter(n_cal==n_cal.plot, n_test == n_test.plot, prop_out<=0.25)
        ## Make plot
        pp <- df |>
        ggplot(aes(x = prop_out, y = Power, color = Method, shape = Method)) +
            geom_line() +
            geom_point() +
            geom_errorbar(aes(ymin = Power - 2*SE, ymax = Power + 2*SE), width = 0.02) +
            geom_hline(yintercept = alpha, linetype = 2) +
            facet_wrap(.~Alternative, nrow=2, labeller="label_value") +
            ylim(0,1) +
            theme_bw(base_size = 15) +
            scale_color_manual(values = colors) +
            scale_shape_manual(values = shapes) +
            labs(#title = paste("Power of Different Methods for Global Testing"),
                 #subtitle = sprintf("Calibration size: %d, Test size: %d", n_cal.plot, n_test.plot),
                 x = "Proportion of Outliers",
                 y = "Power",
                 color = "Method") +
            theme(legend.position = "bottom")       
        ## Save the plot as a PNG file
        filename <- sprintf("figures/power_ncal%d_ntest%d.png", n_cal.plot, n_test.plot)
        ggsave(filename = filename, plot = pp, width = 10, height = 5)
    }

    ##    plot_power_for_n(200, 200)
    plot_power_for_n(500, 200)

    
}

if(plot.power_1) {
    plot_power_1()
}


    
    ## Function to plot and save the power plot for a given alternative value
    plot_power_for_alternative <- function(alt.idx) {
        alternative.plot <- alternative.labels[alt.idx]
        alternative.plot.value <- alternative.values[alt.idx]       
        ## Filter for the specified alternative
        df <- power_results |>
        filter(alternative == alternative.plot)
        ## Make plot
        pp <- df |>
        ggplot(aes(x = prop_out, y = Power, color = Method, shape = Method)) +
            geom_line() +
            geom_point() +
            geom_errorbar(aes(ymin = Power - 2*SE, ymax = Power + 2*SE), width = 0.02) +
            geom_hline(yintercept = alpha, linetype = 2) +
            facet_grid(.~n_cal, labeller="label_both") +
            ylim(0,1) +
            theme_bw(base_size = 15) +
            labs(title = paste("Power of Different Methods for Global Testing"),
                 subtitle = sprintf("Alternative: %s", alternative.plot),
                 x = "Proportion of Outliers",
                 y = "Power",
                 color = "Method") +
            theme(legend.position = "bottom")       
        ## Save the plot as a PNG file
        filename <- paste0("figures/power_", alternative.plot.value, ".png")
        ggsave(filename = filename, plot = pp, width = 8, height = 6)
    }

    ## Make all the plots
    for (alt.idx in seq_along(alternative.labels)) {
        plot_power_for_alternative(alt.idx)
    }

    
    ## Example usage:
    plot_power_for_alternative(1)

    ## Plot the power for different methods and prop_out values
    power_results |>
    filter(n_cal==200, n_test==200) |>
    filter(Method %in% c(methods.plot)) |>
    ggplot(aes(x = prop_out, y = Power, color = Method, shape=Method)) +
        geom_line() +
        geom_point() +
        geom_errorbar(aes(ymin = Power - 2*SE, ymax = Power + 2*SE), width = 0.02) +
        geom_hline(yintercept=alpha, linetype=2) +
        facet_grid(.~alternative) +
        theme_minimal(base_size = 15) +
        labs(title = "Power of Different Methods for Various Proportions of Outliers",
             x = "Proportion of Outliers",
             y = "Power",
             color = "Method")

    
}
