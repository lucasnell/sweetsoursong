
#'
#' Stochastic simulations for supplemental figures.
#'

library(sweetsoursong)
library(tidyverse)
library(patchwork)
library(ggtext)
library(grid)
library(RcppParallel)
setThreadOptions(numThreads = max(defaultNumThreads() - 2L, 1L))

if (file.exists(".Rprofile")) source(".Rprofile")


big_stoch_sims_file <- "_data/big-stoch-sims.rds"
big_open_stoch_sims_file <- "_data/big-stoch-open-sims.rds"


outcome_pal <- c("coexistence" = "#993399",
                 "yeast only" = "#FFCC33",
                 "bacteria only" = "#333399",
                 "extinct" = "gray60")
outcome_shapes <- c("coexistence" = 19,
                    "yeast only" = 19,
                    "bacteria only" = 19,
                    "extinct" = 4)

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "magenta")


# ===========================================================================*
# ===========================================================================*
# functions, objects ----
# ===========================================================================*
# ===========================================================================*


#' These are the bounds of where coexistence starts to occur (calculated
#' in Mathematica), plus the middle of these bounds (when d_b0 = 0.3 and u = 1).
#' We then moved the values that result in exclusion a bit further from the
#' border so that the outcomes don't happen extremely slowly.
d_yp__ = c(`bacteria only` = 0.886308,
           coexistence = 1.17755,
           `yeast only` = 1.46879) +
    c(-1, 0, 1) * 0.1




closed_sims <- function(.u, .d_yp, .n_sigma, .season_surv, .rand_season, .np,
                        .threshold = 1e-6) {

    closed_sim_df <- plant_metacomm_stoch(np = .np, u = .u, d_yp = .d_yp,
                                       n_sigma = .n_sigma,
                                       season_surv = .season_surv,
                                       rand_season = .rand_season)
    season_len <- formals(plant_metacomm_stoch)[["season_len"]]
    n_reps <- formals(plant_metacomm_stoch)[["n_reps"]]

    closed_sim_df |>
        # Take median through time (for each rep and plant) for last season:
        filter(t > max(t) - season_len) |>
        group_by(rep, p) |>
        summarize(Y = median(Y),
                  B = median(B),
                  P = median(P),
                  .groups = "drop") |>
        # Take summed Y and B (across entire landscape) for each rep:
        group_by(rep) |>
        summarize(Y = sum(Y), B = sum(B)) |>
        mutate(outcome = case_when(Y > .threshold & B > .threshold ~ "coexistence",
                                   Y > .threshold & B < .threshold ~ "yeast only",
                                   Y < .threshold & B > .threshold ~ "bacteria only",
                                   TRUE ~ "extinction") |>
                   factor(levels = names(outcome_pal))) |>
        split(~ outcome, drop = FALSE) |>
        imap_dfr(\(x, n) tibble(outcome = n, prob = nrow(x) / n_reps)) |>
        mutate(outcome = factor(outcome, levels = names(outcome_pal))) |>
        mutate(u = .u, d_yp = .d_yp, n_sigma = .n_sigma,
               season_surv = .season_surv,
               rand_season = .rand_season,
               n_plants = .np)
}



open_sims <- function(.u, .d_yp, .n_sigma, .season_surv, .rand_season, .np) {

    open_sim_df <- plant_metacomm_stoch(np = .np, u = .u, d_yp = .d_yp,
                                        n_sigma = .n_sigma,
                                        season_surv = .season_surv,
                                        rand_season = .rand_season,
                                        closed = FALSE)

    open_sim_df |>
        group_by(rep) |>
        summarize(dive = diversity_vector(Y, B, group_size = .np),
                  diss = dissimilarity_vector(Y, B, group_size = .np, TRUE)) |>
        mutate(u = .u, d_yp = .d_yp, n_sigma = .n_sigma,
               season_surv = .season_surv,
               rand_season = .rand_season,
               n_plants = .np)
}




# ===========================================================================*
# ===========================================================================*
# simulations ----
# ===========================================================================*
# ===========================================================================*



if (! file.exists(big_stoch_sims_file)) {

    # Takes ~40 min
    set.seed(1472844374)
    stoch_sim_df <- crossing(.u = 0:10,
                             .d_yp = d_yp__,
                             .n_sigma = c(100, 200, 400),
                             .season_surv = c(0.05, 0.1, 0.2),
                             .rand_season = c(TRUE, FALSE),
                             .np = c(2, 10)) |>
        # Don't do this in parallel bc plant_metacomm_stoch is already
        # doing that
        pmap_dfr(closed_sims)

    write_rds(stoch_sim_df, big_stoch_sims_file)

} else {

    stoch_sim_df <- read_rds(big_stoch_sims_file)

}




if (! file.exists(big_open_stoch_sims_file)) {

    # Takes ~70 min
    set.seed(536842932)
    open_stoch_sim_df <- crossing(.u = 0:10,
                                  .d_yp = d_yp__,
                                  .n_sigma = c(100, 200, 400),
                                  .season_surv = c(0.05, 0.1, 0.2),
                                  .rand_season = c(TRUE, FALSE),
                                  .np = c(2, 10)) |>
        # Don't do this in parallel bc plant_metacomm_stoch is already
        # doing that
        pmap_dfr(open_sims)

    write_rds(open_stoch_sim_df, big_open_stoch_sims_file)

} else {

    open_stoch_sim_df <- read_rds(big_open_stoch_sims_file)

}







# ===========================================================================*
# ===========================================================================*
# plots ----
# ===========================================================================*
# ===========================================================================*


# ------------*
# _ closed ----
# ------------*


#' This creates separate plots by rand_season and n_plants for closed
#' simulations, where I'm going to create group plots by rand_season.
closed_plots <- function(rand_season__, n_plants__) {

    xtitle <- textGrob(paste("'Strength of microbe effects on pollinators",
                             "('*italic(u)*')'") |>
                           (\(x) parse(text = x))(),
                       x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                       gp = gpar(fontsize = 14))

    dzn <- "ABC\nDDD"

    # For labelling within a plot by number of plants:
    letter__ <- letters[which(sort(unique(stoch_sim_df$n_plants)) == n_plants__)]

    map(names(d_yp__), \(n) {
        p <- stoch_sim_df |>
            filter(d_yp == d_yp__[[n]],
                   rand_season == rand_season__,
                   n_plants == n_plants__) |>
            mutate(n_sigma = paste("*n<sub>&sigma;</sub>* =", n_sigma),
                   season_surv = paste("*s* =", season_surv)) |>
            ggplot(aes(u, prob, color = outcome)) +
            ggtitle(paste("Outcome without stochasticity:\n", n)) +
            geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
            geom_point(aes(shape = outcome)) +
            geom_line() +
            facet_grid(season_surv ~ n_sigma) +
            scale_y_continuous("Percent of simulations",
                               limits = c(0, 1),
                               breaks = (0:4 * 25) / 100,
                               labels = c("0%", "", "50%", "", "100%")) +
            scale_x_continuous(breaks = seq(0, 10, 2.5),
                               labels = c("0", "", "5", "", "10")) +
            scale_shape_manual(NULL, values = outcome_shapes) +
            scale_color_manual(NULL, values = outcome_pal) +
            theme(plot.title = element_text(size = 12,
                                            margin = margin(0,0,0,b=6)),
                  axis.title.y = element_text(size = 14),
                  strip.text = element_markdown(),
                  axis.title.x = element_blank())
        if (n != names(d_yp__)[[1]]) {
            p <- p + theme(axis.title.y = element_blank(),
                           axis.text.y = element_blank())
        }
        if (n != tail(names(d_yp__), 1)) {
            p <- p + theme(strip.text.y = element_blank())
        }
        return(p)
    }) |>
        c(list(xtitle, guides = "collect")) |>
        do.call(what = wrap_plots)+
        plot_layout(design = dzn, heights = c(20, 1)) +
        plot_annotation(title = sprintf("(%s)", letter__),
                        subtitle = sprintf("%i plants", n_plants__),
                        theme = theme(plot.title = element_text(size = 20,
                                                                face = "bold",
                                                                hjust = 0),
                                      plot.subtitle = element_markdown(size = 18,
                                                                       hjust = 0.5)))
}


closed_plots_list <- crossing(rand_season__ = unique(stoch_sim_df$rand_season),
                              n_plants__ = unique(stoch_sim_df$n_plants)) |>
    arrange(rand_season__, n_plants__) |>
    pmap(closed_plots)


closed_plots_list |> length()
closed_plots_list[[3]]

closed_plot_combined <- function(i) {
    stopifnot(i %in% 1:2)
    vp1 <- viewport(x = 0.5, y = 1.0, width = 1, height = 0.5,
                    just = c("center", "top"))
    vp2 <- viewport(x = 0.5, y = 0.5, width = 1, height = 0.5,
                    just = c("center", "top"))
    grid.newpage()
    print(closed_plots_list[[(i*2L-1L)]], vp = vp1)
    print(closed_plots_list[[(i*2L)]], vp = vp2)
}


for (i in 1:2) {
    .rs <- sort(unique(stoch_sim_df$rand_season))[[i]]
    fn <- sprintf("_figures/stoch-closed-supp-rand_season%s.pdf", .rs)
    save_plot(fn, closed_plot_combined, 10, 10, fun_args = list(i))
}; rm(i, .rs, fn)




# ------------*
# _ open ----
# ------------*


#' This creates separate plots by rand_season and n_plants for open
#' simulations, where I'm going to create group plots by rand_season.
open_plots <- function(rand_season__, n_plants__) {

    xtitle <- textGrob(paste("'Strength of microbe effects on pollinators",
                             "('*italic(u)*')'") |>
                           (\(x) parse(text = x))(),
                       x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                       gp = gpar(fontsize = 14))

    dzn <- "ABC\nDDD"

    # For labelling within a plot by number of plants:
    letter__ <- letters[which(sort(unique(stoch_sim_df$n_plants)) == n_plants__)]

    map(names(d_yp__), \(n) {
        dd <- open_stoch_sim_df |>
            filter(d_yp == d_yp__[[n]],
                   rand_season == rand_season__,
                   n_plants == n_plants__) |>
            mutate(n_sigma = paste("*n<sub>&sigma;</sub>* =", n_sigma),
                   season_surv = paste("*s* =", season_surv)) |>
            pivot_longer(dive:diss) |>
            mutate(name = factor(name, levels = c("diss", "dive"),
                                 labels = c("dissimilarity", "diversity")))
        dds <- dd |>
            group_by(u, d_yp, n_sigma, season_surv, rand_season, name) |>
            summarize(value = median(value), .groups = "drop")
        # Because maximum metric values vary a lot by d_yp:
        ymax <- c(0.53, 1, 1)[which(names(d_yp__) == n)]
        ymax1 <- round(ymax, 1)
        p <- dd |>
            ggplot(aes(u, value, color = name)) +
            ggtitle(paste("Outcome without stochasticity:\n", n)) +
            geom_point(shape = 1, alpha = 0.1) +
            geom_line(data = dds, linewidth = 1) +
            facet_grid(season_surv ~ n_sigma) +
            scale_y_continuous("Statistic value",
                               limits = c(0, ymax),
                               breaks = (0:4 * ymax1 / 4),
                               labels = paste(c("0.0", "", ymax1/2, "",
                                                sprintf("%.01f", ymax1)))) +
            scale_x_continuous(breaks = seq(0, 10, 2.5),
                               labels = c("0", "", "5", "", "10")) +
            scale_color_viridis_d(NULL, option = "magma", begin = 0.3, end = 0.8) +
            theme(plot.title = element_text(size = 12,
                                            margin = margin(0,0,0,b=6)),
                  axis.title.y = element_text(size = 14),
                  strip.text = element_markdown(),
                  axis.title.x = element_blank())
        if (n != names(d_yp__)[[1]]) {
            p <- p + theme(axis.title.y = element_blank())
        }
        if (n != tail(names(d_yp__), 1)) {
            p <- p + theme(strip.text.y = element_blank())
        }
        return(p)
    }) |>
        c(list(xtitle, guides = "collect")) |>
        do.call(what = wrap_plots)+
        plot_layout(design = dzn, heights = c(20, 1)) +
        plot_annotation(title = sprintf("(%s)", letter__),
                        subtitle = sprintf("%i plants", n_plants__),
                        theme = theme(plot.title = element_text(size = 20,
                                                                face = "bold",
                                                                hjust = 0),
                                      plot.subtitle = element_markdown(size = 18,
                                                                       hjust = 0.5)))
}




open_plots_list <- crossing(rand_season__ = unique(stoch_sim_df$rand_season),
                              n_plants__ = unique(stoch_sim_df$n_plants)) |>
    arrange(rand_season__, n_plants__) |>
    pmap(open_plots)


open_plots_list |> length()
open_plots_list[[3]]


open_plot_combined <- function(i) {
    stopifnot(i %in% 1:2)
    vp1 <- viewport(x = 0.5, y = 1.0, width = 1, height = 0.5,
                    just = c("center", "top"))
    vp2 <- viewport(x = 0.5, y = 0.5, width = 1, height = 0.5,
                    just = c("center", "top"))
    grid.newpage()
    print(open_plots_list[[(i*2L-1L)]], vp = vp1)
    print(open_plots_list[[(i*2L)]], vp = vp2)
}


for (i in 1:2) {
    .rs <- sort(unique(stoch_sim_df$rand_season))[[i]]
    fn <- sprintf("_figures/stoch-open-supp-rand_season%s.pdf", .rs)
    save_plot(fn, open_plot_combined, 10, 10, fun_args = list(i))
}; rm(i, .rs, fn)


