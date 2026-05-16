## -----------------------------------------------------------------------------
## Required results : figA21
## Produces         : Figure A21, Figure A22
## -----------------------------------------------------------------------------

source("plotting_utils_A.R")

## ---- Plot function: lower bound vs proportion of top scores, faceted by dataset
make_plot_lower_bound_data_sel <- function(fig.name, plot.quantile=0.5, save.plot=TRUE) {
    present_methods = c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_auto", "lb_wmw")
    init_settings(present_methods = present_methods)
    summary <- results %>%
        pivot_longer(all_of(present_methods), names_to="method", values_to="lower_bound") %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier, method) %>%
        summarise(LB=quantile(lower_bound, plot.quantile), LB.se=sd(lower_bound)/sqrt(n()))
    df <- summary %>%
        filter(prop_out %in% c(0, 0.2, 0.5)) %>%
        filter(method %in% method.values) %>%
        mutate(Method = factor(method, method.values, method.labels)) %>%
        mutate(Classifier = factor(Classifier, classifier.values, classifier.labels)) %>%
        mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
        mutate(Data = factor(Data, data.values, data.labels))
    df.ref <- results %>%
        filter(prop_out %in% c(0, 0.2, 0.5)) %>%
        group_by(Data, p, Signal, n_train, n_cal, n_test, prop_out, n_out, selection, selected_num, Alpha, Classifier) %>%
        summarise(LB=mean(n_out_sel)) %>%
        mutate(N_out = sprintf("%d outliers", prop_out*n_test)) %>%
        mutate(Data = factor(Data, data.values, data.labels))
    df.ghost <- tibble(prop_out = c(0,0), selected_num=c(0,0), n_test=c(1000,1000), LB=c(0,1)) %>%
        mutate(N_out = sprintf("%d outliers", prop_out*n_test))
    pp <- df %>%
        ggplot(aes(x=selected_num/n_test, y=LB)) +
        geom_point(aes(color=Method, shape=Method, alpha=Method)) +
        geom_line(aes(color=Method, alpha=Method)) +
        geom_line(data=df.ref, linetype=2, aes(x=selected_num/n_test, y=LB), color="black", alpha=1) +
        geom_point(data=df.ghost, aes(x=selected_num/n_test, y=LB), alpha=0) +
        facet_wrap(.~Data, labeller="label_value", scale="free", nrow=2) +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_alpha_manual(values=alpha.scale) +
        xlab("Proportion of top scores in selected test set") +
        ylab("90% lower bound") +
        scale_x_continuous(trans='log10', limits=c(0.01,1)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    if(save.plot) {
        plot.file.1 <- sprintf("figures/%s_data_lower_bound_sel_q%s.pdf",
                               fig.name, plot.quantile)
        ggsave(pp, file=plot.file.1, height=3.5, width=7, units="in")
    } else {
        print(pp)
    }
}

results <- load_data("figA21")

## ---- Run --------------------------------------------------------------------
make_plot_lower_bound_data_sel(fig.name="figA21", plot.quantile=0.5)
make_plot_lower_bound_data_sel(fig.name="figA22", plot.quantile=0.9)
