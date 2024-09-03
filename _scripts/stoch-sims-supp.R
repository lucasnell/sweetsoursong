
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
        summarize(value = mean(value), .groups = "drop")
    p <- dd |>
        ggplot(aes(u, value, color = name)) +
        geom_point(shape = 1, alpha = 0.1) +
        geom_line(data = dds, linewidth = 1) +
        facet_grid(rand_season ~ d_yp) +
        scale_y_continuous("Community metric value (*H* or *BC*)",
                           limits = c(0, 1),
                           breaks = (0:4 * 0.25),
                           labels = paste(c("0.0", "", "0.5", "", "1.0"))) +
        scale_x_continuous(breaks = seq(0, 10, 2.5),
                           labels = c("0", "", "5", "", "10")) +
        scale_color_viridis_d(NULL, option = "viridis", begin = 0.3, end = 0.9) +
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
    fn <- sprintf("_figures/stoch-closed-np=%02i.pdf", np)
    fxn <- \() plot(one_stoch_plot(np, closed__ = TRUE, no_labs = TRUE))
    save_plot(fn, fxn, 4, 2)
}; rm(np, fn, fxn)



#' Open simulations will be shown in supplement so will be combined here:
open_plot <- (one_stoch_plot(2L, closed__ = FALSE, add_title = TRUE,
                             axis.title.y = element_markdown(size = 14),
                             axis.title.x = element_blank()) +
                  geom_text(data = tibble(u = rep(10, 2),
                                          value = c(0.6, 0.92),
                                          d_yp = factor(rep(names(d_yp__)[2],2),
                                                        levels = names(d_yp__)),
                                          rand_season = factor(
                                              c(1, 1), levels = 1:2,
                                              labels = c(
                                                  "non-random<br>starts",
                                                  "random<br>starts")),
                                          name = factor(1:2, labels = c(
                                              "diversity", "dissimilarity"))),
                            aes(label = name), hjust = 1, vjust = c(1, 0),
                            fontface = "bold", size = 10, size.unit = "pt")) +
    one_stoch_plot(10L, closed__ = FALSE, add_title = TRUE,
                   axis.title.y = element_blank(),
                   axis.title.x = element_blank()) +
    textGrob(paste("'Strength of microbe effects on pollinators",
                   "('*italic(u)*')'") |>
                 (\(x) parse(text = x))(),
             x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
             gp = gpar(fontsize = 14)) +
    plot_layout(design = "AB\nCC", heights = c(20, 1)) +
    plot_annotation(tag_levels = list(c("(a)", "(b)", ""))) &
    theme(plot.tag = element_text(size = 16, face = "bold"),
          legend.position = "none")

open_plot

save_plot("_figures/stoch-open.pdf", open_plot, w = 8, h = 5)

