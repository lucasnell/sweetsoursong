
#'
#' This script creates the following files inside `_figures`:
#'
#'  - `invasion-poll-rare-bacteria.pdf`
#'  - `invasion-poll-rare-yeast.pdf`
#'  - `invasion-wins-rare-bacteria.pdf`
#'  - `invasion-wins-rare-yeast.pdf`
#'  - `invasion-occup-rare-bacteria.pdf`
#'  - `invasion-occup-rare-yeast.pdf`
#'  - `supp-u-occup-invasion.pdf`
#'  - `supp-u-density-invasion.pdf`
#'
#' All files are based on stochastic simulations where each species invaded
#' a 100-plant landscape by occupying one plant.
#'
#' Files `invasion-poll-rare-*.pdf` are used inside figure 5 (panels a,d)
#' to show how pollinator visits differs between invader and resident
#' when each species starts out rare.
#'
#' Files `invasion-wins-rare-*.pdf` are used inside figure 5 (panels b,e) to
#' show the proportion of plants where the invader outcompetes the resident
#' at the beginning of the 10th season.
#'
#' Files `invasion-occup-rare-*.pdf` are used inside figure 5 (panels c,f) to
#' show the occupancy of plants (abundance > 0) for the invader at the end of
#' the last (20th) season.
#'
#' File `supp-u-occup-invasion.pdf` is figure S2 and shows the effects of `u`
#' on the number of plants the invading species occupied at the end of the
#' last season.
#'
#' File `supp-u-density-invasion.pdf` is figure S3 and shows the effects of `u`
#' on the landscape-wide abundance of the invading species at the end of the
#' last season.
#'
#'


source("_scripts/02-preamble-stoch.R")


#' Change to `TRUE` if you want to write new plots to files.
#' Otherwise, this script will just create the plot objects for viewing.
.write_plots <- FALSE



# ============================================================================*
# ============================================================================*
# Simulations ----
# ============================================================================*
# ============================================================================*

if (! file.exists(stoch_rds_files$mutual_inv)) {

    # Takes ~26 min with 6 threads
    set.seed(1771251461)
    inv_sim_df <- crossing(.u = 0:10,
                           .d_yp = d_yp__,
                           .q = c(0, 0.25, 0.5, 0.75, 0.95),
                           .rare_sp = c("yeast", "bacteria")) |>
        # Don't do this in parallel bc plant_metacomm_stoch is already
        # doing that
        pmap_dfr(\(.u, .d_yp, .q, .rare_sp) {
            if (.rare_sp == "yeast") {
                .Y0 <- c(rep(0.5, 1L), rep(0,   99L))
                .B0 <- c(rep(0,   1L), rep(0.5, 99L))
            } else {
                .B0 <- c(rep(0.5, 1L), rep(0,   99L))
                .Y0 <- c(rep(0,   1L), rep(0.5, 99L))
            }
            z <- plant_metacomm_stoch(np = 100L, u = .u, d_yp = .d_yp, q = .q,
                                      Y0 = .Y0, B0 = .B0,
                                      no_immig = TRUE,
                                      n_sigma = 100,
                                      season_surv = 0.02,
                                      max_t = 3000,
                                      # only calculate abundance based on the
                                      # last half of the last season:
                                      burnin = 3000 - 75 + 0.1) |>
                select(-P) |>
                pivot_longer(Y:B, names_to = "species", values_to = "density")
            x <- z |>
                filter(t == max(t)) |>
                mutate(occup = density > 0) |>
                group_by(rep, species) |>
                summarize(occup = sum(occup), .groups = "drop")
            y <- z |>
                group_by(rep, species, t) |>
                summarize(density = sum(density), .groups = "drop") |>
                group_by(rep, species) |>
                summarize(density = mean(log1p(density)), .groups = "drop")
            bind_cols(x, density = y$density) |>
                mutate(u = .u, d_yp = .d_yp, q = .q, rare_sp = .rare_sp) |>
                add_factors(.exclude = NULL) |>
                mutate(rare_sp = factor(rare_sp, levels = c("yeast",
                                                              "bacteria"))) |>
                select(u, d_yp, q, rare_sp, everything())
        }, .progress = TRUE)

    write_rds(inv_sim_df, stoch_rds_files$mutual_inv)


} else {

    inv_sim_df <- read_rds(stoch_rds_files$mutual_inv)

}



if (! file.exists(stoch_rds_files$mutual_inv_growth)) {

    # Takes ~3.5 min
    set.seed(32506739)
    inv_growth_df <- crossing(.u = 0:10,
                              .d_yp = d_yp__,
                              .rare_sp = c("yeast", "bacteria")) |>
        pmap_dfr(\(.u, .d_yp, .rare_sp) {
            if (.rare_sp == "yeast") {
                .Y0 <- c(rep(0.5, 1L), rep(0,   99L))
                .B0 <- c(rep(0,   1L), rep(0.5, 99L))
            } else {
                .B0 <- c(rep(0.5, 1L), rep(0,   99L))
                .Y0 <- c(rep(0,   1L), rep(0.5, 99L))
            }
            plant_metacomm_stoch(np = 100L,
                                 u = .u,
                                 Y0 = .Y0, B0 = .B0,
                                 d_yp = .d_yp,
                                 q = 0.5,
                                 n_sigma = 100,
                                 max_t = 3000,
                                 begin_end = TRUE) |>
                mutate(dY = calc_dYdt(Y, P, N = 1 - Y - B, d_yp = .d_yp),
                       dB = calc_dBdt(B, P, N = 1 - Y - B)) |>
                mutate(u = .u, d_yp = .d_yp, rare_sp = .rare_sp)
        }, .progress = TRUE)

    write_rds(inv_growth_df, stoch_rds_files$mutual_inv_growth)


} else {

    inv_growth_df <- read_rds(stoch_rds_files$mutual_inv_growth)

}




# ============================================================================*
# ============================================================================*
# Functions ----
# ============================================================================*
# ============================================================================*

calc_dYdt <- function(Y, P, N, d_yp, L0 = 0.5, m = 0.1) {

    n <- length(Y)

    stopifnot(length(P) == 1 || length(P) == n)
    stopifnot(length(N) == 1 || length(N) == n)
    stopifnot(length(d_yp) == 1 || length(d_yp) == n)
    stopifnot(length(L0) == 1)
    stopifnot(length(m) == 1)

    .delta <- d_yp * (P / (L0 + P)) * Y
    dYdt <- .delta * N - m * Y

    return(dYdt)
}


calc_dBdt <- function(B, P, N,
                      d_b0 = formals(plant_metacomm_stoch)[["d_b0"]],
                      d_bp = formals(plant_metacomm_stoch)[["d_bp"]],
                      L0 = 0.5, m = 0.1) {

    n <- length(B)

    stopifnot(length(P) == 1 || length(P) == n)
    stopifnot(length(N) == 1 || length(N) == n)
    stopifnot(length(d_b0) == 1 || length(d_b0) == n)
    stopifnot(length(d_bp) == 1 || length(d_bp) == n)
    stopifnot(length(L0) == 1)
    stopifnot(length(m) == 1)

    .gamma <- (d_b0 + d_bp * (P / (L0 + P))) * B
    dBdt <- .gamma * N - m * B

    return(dBdt)
}


#'
#' NOTE: The function below calculates P for a specific scenario that was
#' used as the starting conditions in the mutual invasibility simulations.
#' The rare species occupies one plant with density = 0.5
#' (the other with density = 0)
#' The common species occupies 99 plants of the 100 total plants.
#'
calc_P <- function(u, species, rare_sp) {

    #
    # for patch i, equation for P is: n * (1 - B[i])^u / (2 * sum((1 - B[j])^u))
    # Below,
    # Bi = (1 - B[i])^u
    # sumBj = sum((1 - B[j])^u)
    # ... for the scenario described above.
    #
    if (length(species) == 1) species <- rep(species, length(u))
    if (length(rare_sp) == 1) rare_sp <- rep(rare_sp, length(u))

    stopifnot(length(u) == length(species))
    stopifnot(length(u) == length(rare_sp))

    Bi <- numeric(length(u))
    spY <- species == "yeast"
    Bi[spY] <- 1
    Bi[!spY] <- 0.5^(u[!spY])

    sumBj <- numeric(length(u))
    rsY <- rare_sp == "yeast"
    sumBj[rsY] <- 99 * 0.5^(u[rsY]) + 1
    sumBj[!rsY] <- 99 + 0.5^(u[!rsY])

    P <- 100 * Bi / (2 * sumBj)

    return(P)
}

#'
#' Plot occupancy or density after invasion.
#'
invasion_plotter <- function(x_var, y_var, double_facet, .species = NULL) {

    y_var <- match.arg(y_var, c("occup", "density"))

    f_var <- ifelse(x_var == "u", "q", "u")

    if (x_var == "u") {
        x_axis <- scale_x_continuous(paste("Strength of microbe--pollinator",
                                           "effect (*u*)"),
                                     breaks = seq(0, 10, 2.5),
                                     labels = c("0", "", "5", "", "10"))
        dd <- inv_sim_df
    } else {
        x_axis <- scale_x_continuous("Between-season determinism (*q*)",
                                     breaks = as.numeric(levels(inv_sim_df$q)),
                                     labels = c("0", "", "0.5", "", "0.95"))
        dd <- inv_sim_df |> filter(u %in% c(0, 1, 4, 7, 10))
    }

    dd <- dd |>
        mutate(occup = log1p(occup)) |> # density is already transformed
        group_by(u, d_yp, q, rare_sp, species) |>
        summarize(lo = quantile(.data[[y_var]], 0.2),
                  hi = quantile(.data[[y_var]], 0.8),
                  !!y_var := mean(.data[[y_var]]),
                  .groups = "drop")
    y_max <- max(dd[[y_var]])

    if (y_var == "occup") {
        y_axis <- scale_y_continuous("Last-season occupancy",
                                     breaks = log1p(c(0, 3^(0:3))),
                                     labels = c(0, 3^(0:3)),
                                     limits = c(0, y_max))
        y0 <- 1
    } else {
        y_axis <- scale_y_continuous("Last-season abundance",
                                     breaks = log1p(2^(0:3 * 2L - 1L)),
                                     labels = 2^(0:3 * 2L - 1L),
                                     limits = c(0, y_max))
        y0 <- 0.5
    }

    dd <- dd |>
        filter(rare_sp == species) |>
        mutate(increasing = factor(.data[[y_var]] > log1p(y0))) |>
        mutate(!!x_var := as.numeric(paste(.data[[x_var]])))

    if (!is.null(.species)) dd <- filter(dd, species == .species)

    if (double_facet) {
        facet <- facet_grid(vars(.data[[f_var]]), vars(d_yp),
                            labeller = labeller(!!f_var := \(x) lapply(
                                x, \(z) paste0("*", f_var, "* = ", z)),
                                d_yp = label_value))
    } else {
        dd <- dd |>
            filter(.data[[f_var]] == median(as.numeric(levels(.data[[f_var]]))))
        facet <- facet_grid(~ d_yp)
    }

    dd |>
        ggplot(aes(.data[[x_var]], .data[[y_var]], color = species)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_hline(yintercept = log1p(y0), linetype = "22", color = "gray70") +
        geom_point(aes(shape = increasing), size = 2) +
        geom_line(linewidth = 1, show.legend = TRUE) +
        x_axis +
        y_axis +
        facet +
        scale_color_manual(NULL, values = spp_pal) +
        scale_shape_manual(NULL, values = c("TRUE" = 19, "FALSE" = 1),
                           guide = "none") +
        theme(axis.title.x = element_markdown(margin = margin(0,0,0,t=12)),
              strip.text = element_markdown())

}



# ============================================================================*
# ============================================================================*
# Create plots ----
# ============================================================================*
# ============================================================================*


#'
#' Main text figure panels showing pollinator visits by species when
#' each is rare:
#'
inv_poll_ps <- crossing(u = round(seq(0, 10, length.out = 101), 1),
                       species = c("yeast", "bacteria"),
                       rare_sp = c("yeast", "bacteria")) |>
    mutate(P = calc_P(u, species, rare_sp)) |>
    split(~ rare_sp) |>
    map(\(d) {
        d |>
            ggplot(aes(u, P, color = species)) +
            geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
            geom_line(linewidth = 1) +
            scale_color_manual(NULL, values = spp_pal, guide = "none") +
            scale_x_continuous(paste("Strength of microbe--pollinator",
                                     "effect (*u*)"),
                               breaks = seq(0, 10, 2.5),
                               labels = c("0", "", "5", "", "10")) +
            theme(axis.title.y = element_blank(),
                  axis.title.x = element_blank())
    })


if (.write_plots) {
    for (sp in names(inv_poll_ps)) {
        save_plot(sprintf("_figures/invasion-poll-rare-%s.pdf", sp),
                  inv_poll_ps[[sp]], 2, 2)
    }; rm(sp)
} else {
    do.call(wrap_plots, inv_poll_ps) + plot_layout(ncol = 1)
}





#'
#' Proportion where growth for rare species > growth for common species,
#' for each species being rare.
#'
inv_prop_wins_plots <- unique(inv_growth_df$rare_sp) |>
    set_names() |>
    map(\(.rs) {
        inv_growth_df |>
            # beginning of 10th season (round used bc of R rounding issues):
            filter(round(t, 1) == (10 - 1) * 150 + 0.1) |>
            filter(rare_sp == .rs) |>
            group_by(rep, u, d_yp, rare_sp) |>
            summarize(p_inv = ifelse(.rs == "yeast",
                                     mean(dY > dB),
                                     mean(dB > dY)),
                      .groups = "drop") |>
            rename(species = rare_sp) |>
            add_factors(.exclude = "u") |>
            ggplot(aes(u, p_inv, color = species)) +
            geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
            geom_point(shape = 1, size = 2, alpha = 0.1) +
            stat_summary_bin(geom = "line", fun = "mean",
                             breaks = -1:10 + 0.5,
                             linewidth = 1) +
            scale_color_manual(values = spp_pal, guide = "none") +
            scale_y_continuous("Proportion flowers where invader wins",
                               limits = c(0, 1)) +
            scale_x_continuous(paste("Strength of microbe--pollinator",
                                     "effect (*u*)"),
                               breaks = seq(0, 10, 2.5),
                               labels = c("0", "", "5", "", "10")) +
            facet_grid(~ d_yp) +
            theme(axis.title = element_markdown())
    })

if (.write_plots) {
    for (sp in names(inv_prop_wins_plots)) {
        p <- inv_prop_wins_plots[[sp]] + subpanel_theme +
            theme(panel.spacing.x = unit(0.5, "lines"))
        save_plot(sprintf("_figures/invasion-wins-rare-%s.pdf", sp), p, 3.5, 2)
    }; rm(sp, p)
} else {
    do.call(wrap_plots, inv_prop_wins_plots) + plot_layout(ncol = 1)
}



#'
#' Occupancy in final season of the rare species, for each species being rare.
#'
inv_prop_occup_plots <- unique(inv_growth_df$rare_sp) |>
    set_names() |>
    map(\(.rs) {
        invasion_plotter("u", "occup", double_facet = FALSE, .species = .rs)
    })

if (.write_plots) {
    for (sp in names(inv_prop_occup_plots)) {
        p <- inv_prop_occup_plots[[sp]] + subpanel_theme +
            theme(panel.spacing.x = unit(0.5, "lines"))
        save_plot(sprintf("_figures/invasion-occup-rare-%s.pdf", sp), p, 3.5, 2)
    }; rm(sp, p)
} else {
    do.call(wrap_plots, inv_prop_occup_plots) + plot_layout(ncol = 1)
}



#'
#' Supplemental figure of how `u` affects occupancies and abundances from
#' mutual invasion simulations, across multiple values of `q`:
#'
if (.write_plots) {
    for (y in c("occup", "density")) {
        p <- invasion_plotter("u", y, double_facet = TRUE)
        save_plot(sprintf("_figures/supp-u-%s-invasion.pdf", y), p, 4.5, 6)
    }; rm(p, y)
} else {
    map(c("occup", "density"),
        \(y) invasion_plotter("u", y, double_facet = TRUE)) |>
        do.call(what = wrap_plots) +
        plot_layout(nrow = 1)
}



