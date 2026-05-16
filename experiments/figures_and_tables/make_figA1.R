## -----------------------------------------------------------------------------
## Required results : figA1
## Produces         : Figure A1
## -----------------------------------------------------------------------------

library(ggplot2)

results_file <- "../results/figA1/time_results.csv"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

out_df <- read.csv(results_file)
size <- out_df$n
TOTAL_TIME_MS <- as.matrix(out_df[, c("Preparatory_step","Simes","Storey", "Shiraishi","Fisher","WMW")])

df <- data.frame(
    n = rep(size,6),
    Method = factor(rep(c("Preparatory step","Simes","Storey","Shiraishi","Fisher","WMW"),each=5),levels=c("Preparatory step","Simes","Storey","Shiraishi","Fisher","WMW")) ,
    Time = c(TOTAL_TIME_MS[,1],TOTAL_TIME_MS[,2],TOTAL_TIME_MS[,3],TOTAL_TIME_MS[,4],TOTAL_TIME_MS[,5],TOTAL_TIME_MS[,6])
)

source("plotting_utils_A.R")

                                        # Initialize global aesthetic mappings
init_settings()

                                        # Map the labels used in this plot to the method codes in plotting_utils.R
label_to_code <- c(
    "Simes"            = "lb_simes",
    "Storey"           = "lb_storey_simes",
    "Fisher"           = "lb_fisher",
    "WMW"              = "lb_wmw_k2",
    "Shiraishi"        = "lb_lmp",
    "Preparatory step" = "disc_bh"
)

color.by.method <- c(
    lb_simes        = cbPalette[1],
    lb_storey_simes = cbPalette[1],
    lb_fisher       = cbPalette[4],
    lb_wmw_k2       = cbPalette[3],
    lb_wmw_k3       = cbPalette[6],
    lb_lmp          = cbPalette[7],
    lb_auto         = cbPalette[8],
    disc_bh         = cbPalette[9]
)

shape.by.method <- c(
    lb_simes        = 2,
    lb_storey_simes = 6,
    lb_fisher       = 3,
    lb_wmw_k2       = 1,
    lb_wmw_k3       = 0,
    lb_lmp          = 9,
    lb_auto         = 8,
    disc_bh         = 5
)

plot.method.labels <- names(label_to_code)
plot.color.scale <- unname(color.by.method[label_to_code])
plot.shape.scale <- unname(shape.by.method[label_to_code])

df$Method <- factor(df$Method, levels = plot.method.labels)

                                        # Reference trend lines anchored at n = 5000
n_seq <- seq(5000, 25000, length.out = 200)
n0 <- 5000

trend_df <- rbind(
    data.frame(n = n_seq,
               Time = 4 * (n_seq * log(n_seq)) / (n0 * log(n0)),
               Trend = "O(n log n)"),
    data.frame(n = n_seq,
               Time = 3000 * (n_seq / n0)^2,
               Trend = "O(n^2)")
)

p <- ggplot(df, aes(x = n, y = Time, color = Method, shape = Method, group = Method)) +
    geom_line(data = trend_df,
              aes(x = n, y = Time, linetype = Trend, group = Trend),
              color = "grey40", inherit.aes = FALSE, linewidth = 0.5) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    scale_color_manual(values = plot.color.scale, breaks = plot.method.labels, name = "Method") +
    scale_shape_manual(values = plot.shape.scale, breaks = plot.method.labels, name = "Method") +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000),
                  labels = scales::comma,
                  name = "Time (µs)") +
    scale_x_log10(breaks = unique(df$n), labels = scales::comma) +
    scale_linetype_manual(values = c("O(n log n)" = "dashed", "O(n^2)" = "dotted"), name = "Trend") +
    labs(x = "n") +
    theme_bw() +
    theme(legend.position = "right")

ggsave("figures/figA1_time.pdf", plot = p, width = 6, height = 3)
