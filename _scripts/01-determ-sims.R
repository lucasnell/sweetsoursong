
#'
#' This script creates the following files inside `_figures`:
#'
#'   - `determ-u0-d_yp1.4.pdf`
#'   - `determ-u1-d_yp1.4.pdf`
#'   - `determ-u1-d_yp1.9.pdf`
#'   - `determ-u4-d_yp1.9.pdf`
#'
#' They are used inside figure 3 (panels d-g) to show time series for outcomes
#' of competition among yeast and bacteria for the deterministic, two-plant
#' model. These plots differ in their values of `u` (strength of
#' microbe-pollinator effect) and `d_yp` (yeast pollinator-dependent dispersal
#' rate).
#'

source("_scripts/00-preamble.R")

#' Change to `TRUE` if you want to write new plots to files.
#' Otherwise, this script will just create the plot objects for viewing.
.write_plots <- FALSE


#'
#' Simulations to compare dynamics across values of `u` and `d_yp`
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
    if (no_labs) p <- p + subpanel_theme
    if (!is.null(ymax)) p <- p + ylim(0, ymax)
    if (!is.null(.theme)) p <- p + do.call(theme, .theme)
    return(p)
}

u_plots <- u_sims |>
    split(~ u + d_yp, sep = "_", drop = TRUE) |>
    map(one_u_plot, no_labs = TRUE,
        .theme = list(panel.spacing.y = unit(2, "lines")))



if (.write_plots) {
    for (n in names(u_plots)) {
        .u <- strsplit(n, "_")[[1]][[1]]
        .d_yp <- strsplit(n, "_")[[1]][[2]]
        fn <- sprintf("_figures/determ-u%s-d_yp%s.pdf", .u, .d_yp)
        save_plot(fn, u_plots[[n]], 2.25, 2.5)
    }; rm(n, .u, .d_yp, fn)
} else {
    do.call(wrap_plots, u_plots) + plot_layout(nrow = 1)
}

