
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

# data files created here:
data_files <- list(
    # sims of outcomes / abundances:
    main = "_data/big-stoch-sims.rds",
    # sims of mutual invasibility:
    mutual_inv = "_data/big-stoch-mutual-inv-sims.rds",
    # sims of mutual invasiblity - growth rates:
    mutual_inv_growth = "_data/big-stoch-mutual-inv-growth-sims.rds")


outcome_pal <- c("coexistence" = "#008B00",
                 "yeast only" = "#FFCC33",
                 "bacteria only" = "#333399",
                 "extinction" = "gray60")
outcome_shapes <- c("coexistence" = 19,
                    "yeast only" = 19,
                    "bacteria only" = 19,
                    "extinction" = 4)

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "gray60")


# ===========================================================================*
# ===========================================================================*
# data functions, d_yp ----
# ===========================================================================*
# ===========================================================================*


#' "coexistence" is the boundary where coexistence occurs at u=0.
#' The other two are on either side of that boundary and produce competitive
#' exclusion.
d_yp__ = c("bacteria only" = 0.8,
           "coexistence" = 1,
           "yeast only" = 1.2)



add_factors <- function(d, .exclude = NULL) {
    if ("d_yp" %in% colnames(d) && !is.factor(d[["d_yp"]]) &&
        ! "d_yp" %in% .exclude) {
        d <- mutate(d, d_yp = factor(d_yp, levels = unname(d_yp__),
                                     labels = names(d_yp__)))
    }
    if ("q" %in% colnames(d) && !is.factor(d[["q"]]) &&
        ! "q" %in% .exclude) {
        d <- mutate(d, q = factor(q))
    }
    if ("u" %in% colnames(d) && !is.factor(d[["u"]]) &&
        ! "u" %in% .exclude) {
        d <- mutate(d, u = factor(u))
    }
    if ("species" %in% colnames(d) && !is.factor(d[["species"]]) &&
        ! "species" %in% .exclude) {
        if (all(d$species %in% c("Y", "B"))) {
            d <- mutate(d, species = factor(species, levels = c("Y", "B"),
                                            labels = c("yeast", "bacteria")))
        } else if (all(d$species %in% c("yeast", "bacteria"))) {
            d <- mutate(d, species = factor(species,
                                            levels = c("yeast", "bacteria")))
        } else if (all(d$species %in% c("Y", "B", "P"))) {
            d <- mutate(d, species = factor(species, levels = c("Y", "B", "P"),
                                            labels = c("yeast", "bacteria",
                                                       "pollinators")))
        } else if (all(d$species %in% c("yeast", "bacteria", "pollinators"))) {
            d <- mutate(d, species = factor(species,
                                            levels = c("yeast", "bacteria",
                                                       "pollinators")))
        } else {
            stop("strange values in d$species")
        }
    }
    return(d)
}

add_outcome_props <- function(d){

    sc <- c(colnames(d)[!grepl("Y|B|H|rep", colnames(d))], "outcome")
    sf <- paste("~", paste(sc, collapse = " + ")) |> as.formula()

    .threshold <- 0

    d |>
        mutate(outcome = case_when(maxY > .threshold & maxB > .threshold ~
                                       "coexistence",
                                   maxY > .threshold & maxB <= .threshold ~
                                       "yeast only",
                                   maxY <= .threshold & maxB > .threshold ~
                                       "bacteria only",
                                   maxY <= .threshold & maxB <= .threshold ~
                                       "extinction",
                                   TRUE ~ NA_character_)) |>
        mutate(outcome = factor(outcome, levels = names(outcome_pal))) |>
        split(sf, drop = FALSE, sep = "__") |>
        imap_dfr(\(x, n) {
            if (nrow(x) > 0) {
                out <- slice(x, 1) |>
                    select(all_of(sc))
            } else {
                out <- tibble(.rows = 1)
                for (i in 1:length(sc)) out[[sc[i]]] <-
                        str_split(n, "__")[[1]][[i]]
                for (.c in sc) out[[.c]] <-
                        eval(call(paste0("as.", class(x[[.c]])), out[[.c]]))
            }
            out[["prop"]] <- nrow(x) / 100
            return(out)
        })

}



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



# ===========================================================================*
# ===========================================================================*
# even sims ----
# ===========================================================================*
# ===========================================================================*



if (! file.exists(data_files$main)) {

    # Takes ~34 min with 6 threads
    set.seed(1472844374)
    stoch_sim_df <- crossing(.u = 0:10,
                             .d_yp = d_yp__,
                             .season_surv = c(0.01, 0.02, 0.04),
                             .q = c(0, 0.25, 0.5, 0.75, 0.95)) |>
        # Don't do this in parallel bc plant_metacomm_stoch is already
        # doing that
        pmap_dfr(\(.u, .d_yp, .season_surv, .q) {
            plant_metacomm_stoch(np = 100L, u = .u, d_yp = .d_yp,
                                 q = .q,
                                 no_immig = TRUE,
                                 n_sigma = 100,
                                 season_surv = .season_surv,
                                 # only sample last 5 seasons:
                                 burnin = 3000 * 15 / 20,
                                 # summarize by rep:
                                 summarize = "rep") |>
                mutate(u = .u, d_yp = .d_yp, season_surv = .season_surv,
                       q = .q) |>
                select(u, d_yp, season_surv, q, everything())
        }, .progress = TRUE)

    write_rds(stoch_sim_df, data_files$main)

} else {

    stoch_sim_df <- read_rds(data_files$main)

}




# =====================================================*
#           mutual inv. sims ----
# =====================================================*



if (! file.exists(data_files$mutual_inv)) {

    # Takes ~26 min with 6 threads
    set.seed(1771251461)
    inv_sim_df <- crossing(.u = 0:10,
                           .d_yp = d_yp__,
                           .q = sort(unique(stoch_sim_df$q)),
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
                                      season_surv = median(stoch_sim_df$season_surv),
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

    write_rds(inv_sim_df, data_files$mutual_inv)


} else {

    inv_sim_df <- read_rds(data_files$mutual_inv)

}






# ===========================================================================*
# ===========================================================================*
# mutual inv growth sims ----
# ===========================================================================*
# ===========================================================================*


if (! file.exists(data_files$mutual_inv_growth)) {

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
                filter(round(t %% 150, 1) == 0.1) |>
                mutate(dY = calc_dYdt(Y, P, N = 1 - Y - B, d_yp = .d_yp),
                       dB = calc_dBdt(B, P, N = 1 - Y - B)) |>
                mutate(u = .u, d_yp = .d_yp, rare_sp = .rare_sp)
        }, .progress = TRUE)

    write_rds(inv_growth_df, data_files$mutual_inv_growth)


} else {

    inv_growth_df <- read_rds(data_files$mutual_inv_growth)

}





# ===========================================================================*
# ===========================================================================*
# plot functions ----
# ===========================================================================*
# ===========================================================================*



#'
#' Make objects used for two plotting functions below:
#'
make_plot_objects <- function(x_var, facet_q, facet_u,
                              env, facet_fxn = NULL) {

    .f_var <- ifelse(x_var == "u", "q", "u")
    .f_lvls <- ifelse(x_var == "u", list(facet_q), list(facet_u)) |>
        base::`[[`(1)
    if (length(.f_lvls) > 1) {
        # had to define f_var labeller manually bc label_both with a different
        # separator gives weird error:
        f_labeller <- eval(parse(text=paste0("function(x) { lapply(x, ",
                                             "\\(z) paste('*", .f_var,
                                             "* = ', z)) }")))
        if (is.null(facet_fxn)) facet_fxn <- facet_grid
        .facet <- facet_fxn(vars(.data[[.f_var]]), vars(d_yp),
                             labeller = labeller(!!.f_var := f_labeller,
                                                 d_yp = label_value))
    } else {
        if (is.null(facet_fxn)) facet_fxn <- facet_wrap
        .facet <- facet_fxn( ~ d_yp, nrow = 1)
    }
    if (x_var == "u") {
        .x_scale <- scale_x_continuous(paste("Strength of microbe--pollinator",
                                            "effect (*u*)"),
                                      breaks = seq(0, 10, 2.5),
                                      labels = c("0", "", "5", "", "10"))
    } else {
        .x_scale <- scale_x_continuous("Between-season determinism (*q*)",
                                      breaks = c(0:3 * 0.25, 0.95),
                                      labels = c("0", "", "0.5", "", "0.95"))
    }
    env[["x_scale"]] <- .x_scale
    env[["facet"]] <- .facet
    env[["f_lvls"]] <- .f_lvls
    env[["f_var"]] <- .f_var
    env[["f_lvls"]] <- .f_lvls

    invisible(NULL)
}



outcome_plotter <- function(x_var,
                            facet_q = c(0, 0.25, 0.5, 0.75, 0.95),
                            facet_u = c(0, 1, 4, 7, 10),
                            drop_ext = TRUE,
                            add_title = FALSE) {

    x_var <- match.arg(x_var, c("q", "u"))

    # makes objects: x_scale, facet, f_lvls, f_var, f_lvls
    make_plot_objects(x_var, facet_q, facet_u,
                      env = sys.frame(sys.nframe()))

    dd <- stoch_sim_df |>
        filter(season_surv == median(season_surv),
               .data[[f_var]] %in% f_lvls) |>
        select(q, u, d_yp, maxY, maxB) |>
        add_outcome_props() |>
        add_factors(.exclude = x_var)

    if (drop_ext && sum(dd$prop[dd$outcome == "extinction"]) == 0) {
            dd <- dd |>
                filter(outcome != "extinction") |>
                mutate(outcome = fct_drop(outcome, "extinction"))
    }

    p <- dd |>
        ggplot(aes(.data[[x_var]], prop, color = outcome)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
        geom_point(aes(shape = outcome), show.legend = TRUE) +
        geom_line(show.legend = TRUE) +
        facet +
        scale_y_continuous("Percent of simulations",
                           limits = c(0, 1),
                           breaks = (0:4 * 25) / 100,
                           labels = scales::label_percent()) +
        x_scale +
        scale_shape_manual(NULL, values = outcome_shapes, drop = FALSE) +
        scale_color_manual(NULL, values = outcome_pal, drop = FALSE) +
        theme(axis.title.x = element_markdown(margin = margin(0,0,0,t=12)),
              strip.text.y = element_markdown(),
              plot.title = element_markdown())

    if (add_title) {
        .title <- sprintf("*s* = %s", median(stoch_sim_df$season_surv))
        if (length(f_lvls) == 1) {
            .title <- paste0(.title, sprintf(", *%s* = %s", f_var, f_lvls))
        }
        p <- p + ggtitle(.title)
    }

    return(p)
}





abundance_plotter <- function(x_var,
                              facet_q = c(0, 0.25, 0.5, 0.75, 0.95),
                              facet_u = c(0, 1, 4, 7, 10),
                              log_trans = TRUE,
                              spp__ = NULL,
                              free_y = FALSE,
                              add_title = FALSE) {

    x_var <- match.arg(x_var, c("q", "u"))

    yvar <- ifelse(log_trans, "meanlogY", "meanY")
    bvar <- ifelse(log_trans, "meanlogB", "meanB")

    if (free_y) {
        .facet_fxn <- \(...) facet_wrap(scales = "free_y",
                                        ncol = length(unique(stoch_sim_df$d_yp)),
                                        ...)
        .ylimits <- NULL
        .ybreaks <- waiver()
        .ylabels <- waiver()
        .hline <- list()
    } else {
        .facet_fxn <- NULL
        .ylimits <- ifelse(log_trans, list(log1p(c(0, 100))), list(c(0, 100))) |>
            base::`[[`(1)
        .ybreaks <- ifelse(log_trans, list(log1p(4^(0:3))), list(0:4 * 25)) |>
            base::`[[`(1)
        .ylabels <- ifelse(log_trans, list(c(4^(0:3))),
                           list(c("0", "", "50", "", "100")))[[1]]
        h_line <- geom_hline(yintercept = 0, linetype = 1, color = "gray80")
    }

    # makes objects: x_scale, facet, f_lvls, f_var, f_lvls
    make_plot_objects(x_var, facet_q, facet_u,
                      facet_fxn = .facet_fxn,
                      env = sys.frame(sys.nframe()))

    dd <- stoch_sim_df |>
        filter(season_surv == median(season_surv),
               .data[[f_var]] %in% f_lvls) |>
        select(q, u, d_yp, {{ yvar }}, {{ bvar }}) |>
        pivot_longer({{ yvar }}:{{ bvar }}, names_to = "species") |>
        mutate(species = str_sub(species, nchar(species), nchar(species))) |>
        add_factors(.exclude = x_var)
    dds <- dd |>
        group_by(q, u, d_yp, species) |>
        summarize(lo = min(value),
                  hi = max(value),
                  lo10 = quantile(value, 0.1),
                  hi10 = quantile(value, 0.9),
                  value = mean(value), .groups = "drop")
    if (!is.null(spp__)) {
        dd <- dd |> filter(species == spp__)
        dds <- dds |> filter(species == spp__)
    }
    p <- dd |>
        ggplot(aes(.data[[x_var]], value, color = species)) +
        facet +
        h_line +
        geom_point(aes(color = species), alpha = 0.1, shape = 1) +
        geom_line(data = dds, linewidth = 1) +
        ## If you want to switch back to using ribbons
        ## (also use `dds` instead of `dd`):
        # geom_ribbon(aes(fill = species, ymin = lo, ymax = hi),
        #             alpha = 0.25, color = NA) +
        # geom_ribbon(aes(fill = species, ymin = lo10, ymax = hi10),
        #             alpha = 0.25, color = NA) +
        # geom_line(linewidth = 1) +
        scale_y_continuous("Microbial abundance", limits = .ylimits,
                           breaks = .ybreaks, labels = .ylabels) +
        x_scale +
        scale_shape_manual(NULL, values = c(0, 2), drop = FALSE) +
        scale_linetype_manual(NULL, values = c(1, 1), drop = FALSE) +
        scale_color_manual(NULL, values = spp_pal, drop = FALSE,
                           aesthetics = c("color", "fill")) +
        theme(axis.title.x = element_markdown(margin = margin(0,0,0,t=12)),
              strip.text = element_markdown(),
              plot.title = element_markdown())

    if (add_title) {
        .title <- sprintf("*s* = %s", median(stoch_sim_df$season_surv))
        if (length(f_lvls) == 1) {
            .title <- paste0(.title, sprintf(", *%s* = %s", f_var, f_lvls))
        }
        if (!is.null(spp__)) {
            .title <- paste0(.title, ", ", spp__)
        }
        p <- p + ggtitle(.title)
    }

    return(p)
}





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






# =====================================================*
#           even t0 plots ----
# =====================================================*

# Plots for when species start evenly distributed.

#'
#' Main text figure panels with `q = 0.5`:
#'
main_plots <- list(outcome_plotter("u", facet_q = 0.5, facet_u = 4),
                   abundance_plotter("u", facet_q = 0.5, facet_u = 4)) |>
    map(\(p) {
        p +
            theme(plot.title = element_blank(),
                  legend.position = "none",
                  panel.spacing.x = unit(1, "lines"),
                  axis.title.y = element_blank(),
                  axis.title.x = element_blank(),
                  strip.text = element_blank())
    }) |>
    do.call(what = wrap_plots) +
    plot_layout(ncol = 1)

save_plot("_figures/u-stochastic-main.pdf", main_plots, 4, 3.4)




#'
#' Supplemental figure of outcomes and abundances from evenly distributed
#' starting abundances, with all `q x u` combos:
#'
supp_u_out_abund_p <-
    (outcome_plotter("u", drop_ext = FALSE) +
         theme(plot.title = element_blank())) +
    (abundance_plotter("u") + theme(plot.title = element_blank())) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))

save_plot("_figures/supp-u-out-abund.pdf", supp_u_out_abund_p, 9, 6)




# =====================================================*
#           mutual inv. plots ----
# =====================================================*


# Plots for when one species is rare.

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

do.call(wrap_plots, inv_poll_ps) +
    plot_layout(ncol = 1)


for (sp in names(inv_poll_ps)) {
    save_plot(sprintf("_figures/invasion-poll-rare-%s.pdf", sp),
              inv_poll_ps[[sp]], 2, 2)
}; rm(sp)



inv_prop_win_plots <- unique(inv_growth_df$rare_sp) |>
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

do.call(wrap_plots, inv_prop_win_plots) +
    plot_layout(ncol = 1)




for (sp in names(inv_prop_win_plots)) {
    ps <- list(inv_prop_win_plots[[sp]],
               invasion_plotter("u", "occup", double_facet = FALSE, .species = sp)) |>
        map(\(p) {
            p + theme(axis.title.y = element_blank(),
                      axis.title.x = element_blank(),
                      legend.position = "none",
                      strip.text = element_blank(),
                      panel.spacing.x = unit(0.5, "lines"))
        })
    save_plot(sprintf("_figures/invasion-wins-rare-%s.pdf", sp), ps[[1]], 3.5, 2)
    save_plot(sprintf("_figures/invasion-occup-rare-%s.pdf", sp), ps[[2]], 3.5, 2)
}; rm(sp, ps)







#'
#' Supplemental figure of occupancies and abundances from mutual invasion
#' simulations, with all `q x u` combos:
#'
for (x in c("u", "q")) {
    for (y in c("occup", "density")) {
        p <- invasion_plotter(x, y, double_facet = TRUE)
        save_plot(sprintf("_figures/supp-%s-%s-invasion.pdf", x, y), p, 4.5, 6)
    }
}; rm(p, x, y)





