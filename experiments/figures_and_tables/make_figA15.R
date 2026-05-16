## -----------------------------------------------------------------------------
## Required results : figA15
## Produces         : Figure A15, Figure A16, Figure A17, Figure A18, Figure A19, Figure A20
## -----------------------------------------------------------------------------

source("utils_plotting_A.R")

## ---- Plot function: LB median / 90th quantile / Power, faceted Key x Classifier
make_plot_data_4 <- function(fig.name, plot.data, plot.n_train, plot.n_cal, plot.n_test, save.plot=TRUE) {
    init_settings(idx.exclude=c(6))
    key.values <- c("LB.50", "LB.90", "Power")
    key.labels <- c("LB (median)", "LB (90th q.)", "Power (global)")
    summary <- results %>%
        pivot_longer(c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw", "lb_higher_k2"),
                     names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, Alpha, Classifier, method) %>%
        summarise(Value_LB.50=median(lower_bound),
                  Value_LB.90=quantile(lower_bound, 0.9),
                  Value_Power=mean(lower_bound>0),
                  SE_LB.50=sd(lower_bound)/sqrt(n()),
                  SE_LB.90=sd(lower_bound)/sqrt(n()),
                  SE_Power=sd(lower_bound>0)/sqrt(n()),
                  N=n()) %>%
        pivot_longer(c("Value_LB.50", "Value_LB.90", "Value_Power",
                       "SE_LB.50", "SE_LB.90", "SE_Power"),
                     names_to = c(".value", "Key"), names_sep = "_")
    df <- summary %>%
        filter(Data==plot.data, n_train==plot.n_train, n_cal==plot.n_cal, n_test==plot.n_test) %>%
        filter(method %in% method.values) %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Classifier = factor(Classifier, classifier.values, classifier.labels)) %>%
        mutate(Key = factor(Key, key.values, key.labels))
    x.max <- max(df$n_out)
    df.ref <- tibble(Key=c("LB.50", "LB.50", "LB.90", "LB.90", "Power", "Power"),
                     n_out=c(0,x.max,0,x.max,0,x.max),
                     Value=c(0,x.max,0,x.max,0.1,0.1),
                     method="lb_auto", Target="Closed testing") %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Key = factor(Key, key.values, key.labels))
    df.range <- tibble(Key=c("LB.50", "LB.50", "LB.90", "LB.90", "Power", "Power"),
                       n_out=c(0,x.max,0,x.max,0,x.max),
                       Value=c(0,x.max,0,x.max,0,1),
                       method="lb_auto", Target="Closed testing") %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Key = factor(Key, key.values, key.labels))
    pp <- df %>%
        ggplot(aes(x=n_out, y=Value, color=Method, shape=Method, alpha=Method)) +
        geom_point() +
        geom_line() +
        geom_line(data=df.ref, aes(x=n_out, y=Value), linetype=2, color="black", alpha=1) +
        geom_point(data=df.range, aes(x=n_out, y=Value), color="black", alpha=0) +
        facet_grid(Key~Classifier,
                   labeller=labeller(Classifier = label_both, Key = label_value),
                   scales="free") +
        xlab("True number of outliers") +
        ylab("") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_alpha_manual(values=alpha.scale) +
        xlim(0,x.max) +
        theme_bw() +
        guides(linetype = "none",
               color=guide_legend(title="Local tests"),
               shape=guide_legend(title="Local tests"),
               alpha=guide_legend(title="Local tests"))
    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_setup4_%s_%d_%d_%d.pdf",
                               fig.name, plot.data, plot.n_train, plot.n_cal, plot.n_test)
        ggsave(pp, file=plot.file.1, height=3.5, width=7, units="in")
    } else {
        print(pp)
    }
}

results <- load_data("figA15")

## ---- Run --------------------------------------------------------------------
## One figure per real-data benchmark (A15-A20).
fig.map <- list("figA15" = "creditcard",
                "figA16" = "pendigits",
                "figA17" = "cover",
                "figA18" = "shuttle",
                "figA19" = "mammography",
                "figA20" = "aloi")

plot.n_train <- 1000
plot.n_cal <- 200
plot.n_test <- 100

for(fig.name in names(fig.map)) {
    make_plot_data_4(fig.name=fig.name, plot.data=fig.map[[fig.name]],
                     plot.n_train=plot.n_train, plot.n_cal=plot.n_cal,
                     plot.n_test=plot.n_test, save.plot=TRUE)
}
