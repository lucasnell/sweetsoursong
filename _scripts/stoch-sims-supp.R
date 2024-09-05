
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


outcome_pal <- c("coexistence" = "#008B00",
                 "yeast only" = "#FFCC33",
                 "bacteria only" = "#333399",
                 "extinct" = "gray60")
outcome_shapes <- c("coexistence" = 19,
                    "yeast only" = 19,
                    "bacteria only" = 19,
                    "extinct" = 4)

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "gray60")


# ===========================================================================*
# ===========================================================================*
# functions, objects ----
# ===========================================================================*
# ===========================================================================*


#' These are the bounds of where coexistence starts to occur (calculated
#' in Mathematica), plus the middle of these bounds (when d_b0 = 0.3 and u = 1).
#' We then moved the values that result in exclusion a bit further from the
#' border so that the outcomes don't happen extremely slowly.
d_yp__ = c("bacteria only" = 0.886308,
           "coexistence" = 1.17755,
           "yeast only" = 1.46879) +
    c(-1, 0, 1) * 0.1




# Summarize by different outcomes
summ_outcomes <- function(sim_df, .threshold = 1e-6) {

    season_len <- formals(plant_metacomm_stoch)[["season_len"]]
    n_reps <- formals(plant_metacomm_stoch)[["n_reps"]]

    sim_df |>
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
        mutate(outcome = factor(outcome, levels = names(outcome_pal)))
}


# Summarize by diversity and dissimilarity metrics
summ_metrics <- function(sim_df) {

    .np <- length(levels(sim_df$p))

    sim_df |>
        group_by(rep) |>
        summarize(dive = diversity_vector(Y, B, group_size = .np),
                  diss = dissimilarity_vector(Y, B, group_size = .np, TRUE))
}




# ===========================================================================*
# ===========================================================================*
# simulations ----
# ===========================================================================*
# ===========================================================================*



if (! file.exists(big_stoch_sims_file)) {

    # Takes ~15 min
    set.seed(1472844374)
    stoch_sim_df <- crossing(.u = 0:10,
                             .d_yp = d_yp__,
                             .rand_season = c(TRUE, FALSE),
                             .np = c(2L, 10L),
                             .closed = c(TRUE, FALSE)) |>
        # Don't do this in parallel bc plant_metacomm_stoch is already
        # doing that
        pmap_dfr(\(.u, .d_yp, .rand_season, .np, .closed) {
            z <- plant_metacomm_stoch(np = .np, u = .u, d_yp = .d_yp,
                                      rand_season = .rand_season,
                                      closed = .closed) |>
                filter(t > max(t) / 2)
            tibble(u = .u, d_yp = .d_yp, rand_season = .rand_season,
                   n_plants = .np, closed = .closed,
                   outs = list(summ_outcomes(z)),
                   mets = list(summ_metrics(z)))
        }, .progress = TRUE)

    write_rds(stoch_sim_df, big_stoch_sims_file)

} else {

    stoch_sim_df <- read_rds(big_stoch_sims_file)

}







# ===========================================================================*
# ===========================================================================*
# plots ----
# ===========================================================================*
# ===========================================================================*




one_stoch_plot <- function(n_plants__,
                           closed__,
                           no_labs = FALSE,
                           add_title = FALSE,
                           ribbon = TRUE,
                           ...) {
    dd <- stoch_sim_df |>
        filter(n_plants == n_plants__,
               closed == closed__) |>
        select(rand_season, u, d_yp, mets) |>
        unnest(mets) |>
        mutate(d_yp = factor(d_yp, levels = unname(d_yp__),
                             labels = names(d_yp__)))|>
        pivot_longer(dive:diss) |>
        mutate(name = factor(name, levels = c("dive", "diss"),
                             labels = c("diversity", "dissimilarity")),
               rand_season = factor(rand_season, levels = c(FALSE, TRUE),
                                    labels = c("non-random<br>starts",
                                               "random<br>starts")))
    dds <- dd |>
        group_by(rand_season, u, d_yp, name) |>
        summarize(lo = min(value),  # quantile(value, 0.05),
                  hi = max(value),  # quantile(value, 0.95),
                  value = mean(value), .groups = "drop")
    if (ribbon) {
        p <- dds |>
            ggplot(aes(u, value, color = name)) +
            geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
            geom_ribbon(aes(ymin = lo, ymax = hi, fill = name),
                        color = NA, alpha = 0.25) +
            geom_line(linewidth = 1)
    } else {
        p <- dd |>
            ggplot(aes(u, value, color = name)) +
            geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
            geom_point(shape = 1, alpha = 0.1) +
            geom_line(data = dds, linewidth = 1)
    }
    p <- p +
        facet_grid(rand_season ~ d_yp) +
        scale_y_continuous("Community metric value (*H* or *BC*)",
                           limits = c(0, 1),
                           breaks = (0:4 * 0.25),
                           labels = paste(c("0.0", "", "0.5", "", "1.0"))) +
        scale_x_continuous(breaks = seq(0, 10, 2.5),
                           labels = c("0", "", "5", "", "10")) +
        scale_color_viridis_d(NULL, option = "plasma", begin = 0.1, end = 0.5,
                              aesthetics = c("color", "fill")) +
        theme(strip.text = element_markdown(),
              axis.title = element_markdown(),
              plot.title = element_text(size = 14,
                                        margin = margin(0,0,0,b=9)),
              ...)
    if (no_labs) p <- p + theme(axis.title = element_blank(),
                                strip.text = element_blank(),
                                legend.position = "none")
    if (add_title) p <- p + ggtitle(sprintf("%i plants", n_plants__))
    return(p)
}




#' Closed simulations will be shown in main text and final figure will be
#' edited in Illustrator, so these need to be simpler:
for (np in unique(stoch_sim_df$n_plants)) {
    for (cl in c(TRUE, FALSE)) {
        fn <- sprintf("_figures/stoch-%s-np=%02i.pdf",
                      ifelse(cl, "closed", "open"), np)
        fxn <- \() plot(one_stoch_plot(np, closed__ = cl, no_labs = TRUE))
        save_plot(fn, fxn, 4, 2)
    }
}; rm(np, cl, fn, fxn)






set.seed(1701361777)
ex_nonrand_sims <- plant_metacomm_stoch(np = 2, u = 0, d_yp = d_yp__[[3]],
                                  rand_season = FALSE, closed = TRUE,
                                  max_t = 300)
set.seed(381552468)
ex_rand_sims <- plant_metacomm_stoch(np = 2, u = 0, d_yp = d_yp__[[3]],
                     rand_season = TRUE, closed = TRUE,
                     max_t = 300)

ex_nonrand_sims |>
    filter(t == 150 | t == 150.1) |>
    select(rep, p, t, Y, B) |>
    pivot_longer(Y:B, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B"),
                         labels = c("yeast", "bacteria"))) |>
    pivot_wider(names_from = t, values_from = density, names_prefix = "t") |>
    ggplot(aes(t150, t150.1)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
    geom_vline(xintercept = 0, linetype = 1, color = "gray80") +
    geom_abline(intercept = 0, slope = 1/10, linetype = 2, color = "gray80") +
    geom_point(aes(color = type), size = 2) +
    scale_color_manual(values = spp_pal, guide = "none") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
ex_rand_sims |>
    filter(t == 150 | t == 150.1) |>
    select(rep, p, t, Y, B) |>
    pivot_longer(Y:B, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B"),
                         labels = c("yeast", "bacteria"))) |>
    pivot_wider(names_from = t, values_from = density, names_prefix = "t") |>
    ggplot(aes(t150, t150.1)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
    geom_vline(xintercept = 0, linetype = 1, color = "gray80") +
    geom_abline(intercept = 0, slope = 1/10, linetype = 2, color = "gray80") +
    geom_point(aes(color = type), size = 2) +
    scale_color_manual(values = spp_pal, guide = "none") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

dt <- 2

ex_nonrand_p <- ex_nonrand_sims |>
    filter(p == "patch 2") |>
    filter(t >= (150 - dt), t <= (150 + dt)) |>
    select(-P) |>
    pivot_longer(Y:B, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B"),
                         labels = c("yeast", "bacteria")),
           id = interaction(type, p, rep, drop = TRUE)) |>
    ggplot(aes(t, density, group = id)) +
    geom_hline(yintercept = 0, linewidth = 1, color = "gray80") +
    geom_line(aes(color = type), linewidth = 1, alpha = 0.1) +
    xlab("Time (days)") +
    scale_color_manual(NULL, values = spp_pal, guide = "none") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
ex_rand_p <- ex_rand_sims |>
    filter(p == "patch 2") |>
    filter(t >= (150 - dt), t <= (150 + dt)) |>
    select(-P) |>
    pivot_longer(Y:B, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B"),
                         labels = c("yeast", "bacteria")),
           id = interaction(type, p, rep, drop = TRUE)) |>
    ggplot(aes(t, density, group = id)) +
    geom_hline(yintercept = 0, linewidth = 1, color = "gray80") +
    geom_line(aes(color = type), linewidth = 1, alpha = 0.1) +
    xlab("Time (days)") +
    scale_color_manual(NULL, values = spp_pal, guide = "none") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())



# save_plot("_figures/stoch-season-nonrand.pdf", ex_nonrand_p, w = 2.5, h = 2)
# save_plot("_figures/stoch-season-rand.pdf", ex_rand_p, w = 2.5, h = 2)
