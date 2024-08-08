
library(sweetsoursong)
library(tidyverse)


outcome_pal <- c("coexist" = "#993399",
                 "yeast only" = "#FFCC33",
                 "bacteria only" = "#333399")

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "gray60")




#'
#' Simulations to compare dynamics across values of u and yeast dispersal
#' where one plant has lots of yeast and little bacteria, and another
#' plant is the opposite.
#'
u_sims <- tibble(.u = c(0, 1, 1, 4), .d_yp = c(rep(1.4, 2), rep(1.9, 2))) |>
    pmap_dfr(\(.u, .d_yp) {
        high0 <- 0.5
        low0 <- 0.02
        plant_metacomm(np = 2, u = .u, d_yp = .d_yp,
                       Y0 = c(high0, low0), B0 = c(low0, high0), max_t = 250) |>
            mutate(u = .u, d_yp = .d_yp) |>
            select(u, d_yp, everything())
    }) |>
    pivot_longer(Y:P, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B", "P"),
                         labels = c("yeast", "bacteria", "pollinators")),
           u = factor(u),
           d_yp = factor(d_yp))



one_u_plot <- function(x, ymax = 1, no_labs = FALSE, .theme = NULL) {
    p <- x |>
        filter(t %% 1 == 0) |>
        ggplot(aes(t, density)) +
        geom_hline(yintercept = 0, linewidth = 1, color = "gray80") +
        geom_line(aes(color = type, linetype = type), linewidth = 1) +
        facet_grid(p ~ .) +
        xlab("Time (days)") +
        scale_color_manual(NULL, values = spp_pal) +
        scale_linetype_manual(NULL, values = c("solid", "solid", "24"))
    if (no_labs) p <- p +
            theme(strip.text = element_blank(),
                  strip.background = element_blank(),
                  legend.position = "none",
                  axis.title = element_blank())
    if (!is.null(ymax)) p <- p + ylim(0, ymax)
    if (!is.null(.theme)) p <- p + do.call(theme, .theme)
    return(p)
}

u_plots <- u_sims |>
    split(~ u + d_yp, sep = "_", drop = TRUE) |>
    map(one_u_plot, no_labs = TRUE,
        .theme = list(panel.spacing.y = unit(2, "lines")))


for (n in names(u_plots)) {
    .u <- strsplit(n, "_")[[1]][[1]]
    .d_yp <- strsplit(n, "_")[[1]][[2]]
    fn <- sprintf("_figures/determ-u%s-d_yp%s.pdf", .u, .d_yp)
    save_plot(fn, u_plots[[n]], 2.25, 2.5)
}; rm(n, .u, .d_yp, fn)

