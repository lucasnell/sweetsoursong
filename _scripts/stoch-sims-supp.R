
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

# sims of outcomes / abundances:
main_stoch_sims_file <- "_data/big-stoch-sims.rds"
# sims of gain / loss:
gain_loss_stoch_sims_file <- "_data/big-stoch-gain-loss-sims.rds"


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
# functions, objects ----
# ===========================================================================*
# ===========================================================================*


#' These are the bounds of where coexistence starts to occur (calculated
#' in Mathematica), plus the middle of these bounds (when d_b0 = 0.3 and u = 1).
#' We then moved the values that result in exclusion a bit further from the
#' border so that the outcomes don't happen extremely slowly.
d_yp__ = c("bacteria only" = 0.886308 - 0.1,
           "coexistence" = 1.17755,
           "yeast only" = 1.46879 + 0.1)





# ===========================================================================*
# ===========================================================================*
# simulations, outcomes/abundances ----
# ===========================================================================*
# ===========================================================================*



if (! file.exists(main_stoch_sims_file)) {

    # Takes ~34 min with 6 threads
    set.seed(1472844374)
    stoch_sim_df <- crossing(.u = 0:10,
                             .d_yp = d_yp__,
                             .season_surv = c(0.01, 0.02, 0.04),
                             .q = c(0, 0.25, 0.5, 0.75, 0.95, 1)) |>
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

    write_rds(stoch_sim_df, main_stoch_sims_file)

} else {

    stoch_sim_df <- read_rds(main_stoch_sims_file)

}

#' q = 1 causes problems bc it causes any plant-level extinctions
#' to be permanent
stoch_sim_df <- stoch_sim_df |>
    filter(q < 1)




# ===========================================================================*
# ===========================================================================*
# plots ----
# ===========================================================================*
# ===========================================================================*

add_factors <- function(d, .exclude = "u") {
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

add_outcome_props <- function(d, .threshold = 0){

    sc <- c(colnames(d)[!grepl("Y|B|H|rep", colnames(d))], "outcome")
    sf <- paste("~", paste(sc, collapse = " + ")) |> as.formula()

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


# used for a bunch of functions below
# had to define these bc label_both with a different separator gives weird error
u_labeller <- function(labels) {
    lapply(labels, \(z) paste("*u* = ", z))
}
q_labeller <- function(labels) {
    lapply(labels, \(z) paste("*q* = ", z))
}


# =====================================================*
#           outcomes ----
# =====================================================*



outcome_plotter <- function(season_surv__, q__ = NULL, .threshold = 0) {

    if (is.null(q__)) {
        p <- stoch_sim_df |>
            filter(season_surv == season_surv__) |>
            select(q, u, d_yp, maxY, maxB) |>
            add_outcome_props(.threshold) |>
            add_factors() |>
            ggplot(aes(u, prop, color = outcome)) +
            ggtitle(sprintf("*s* = %s", season_surv__)) +
            facet_grid(q ~ d_yp,
                       labeller = labeller(q = q_labeller, d_yp = label_value))
    } else {
        p <- stoch_sim_df |>
            filter(season_surv == season_surv__,
                   q == q__) |>
            add_factors() |>
            select(u, d_yp, maxY, maxB) |>
            add_outcome_props(.threshold) |>
            ggplot(aes(u, prop, color = outcome)) +
            ggtitle(sprintf("*s* = %s, *q* = %s", season_surv__, q__)) +
            facet_wrap( ~ d_yp, nrow = 1)
    }
    p <- p +
        geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
        geom_point(aes(shape = outcome)) +
        geom_line() +
        scale_y_continuous("Percent of simulations",
                           limits = c(0, 1),
                           breaks = (0:4 * 25) / 100,
                           labels = c("0%", "", "50%", "", "100%")) +
        scale_x_continuous("Strength of microbe--pollinator effect (*u*)",
                           breaks = seq(0, 10, 2.5),
                           labels = c("0", "", "5", "", "10")) +
        scale_shape_manual(NULL, values = outcome_shapes, drop = FALSE) +
        scale_color_manual(NULL, values = outcome_pal, drop = FALSE) +
        theme(axis.title.x = element_markdown(margin = margin(0,0,0,t=12)),
              strip.text.y = element_markdown(),
              plot.title = element_markdown())
    if (season_surv__ != sort(unique(stoch_sim_df$season_surv))[2]) {
        p <- p + theme(axis.title.x = element_blank())
    }
    return(p)
}



outcome_plotter2 <- function(season_surv__,
                             u__ = c(0, 1, 4, 7, 10),
                             .threshold = 0) {

    if (length(u__) > 1) {
        p <- stoch_sim_df |>
            filter(season_surv == season_surv__) |>
            filter(u %in% u__) |>
            select(q, u, d_yp, maxY, maxB) |>
            add_outcome_props(.threshold) |>
            add_factors(.exclude = NULL) |>
            ggplot(aes(q, prop, color = outcome)) +
            ggtitle(sprintf("*s* = %s", season_surv__)) +
            facet_grid(u ~ d_yp,
                       labeller = labeller(u = u_labeller, d_yp = label_value))
    } else {
        p <- stoch_sim_df |>
            filter(season_surv == season_surv__,
                   u == u__) |>
            select(q, d_yp, maxY, maxB) |>
            add_outcome_props(.threshold) |>
            # add_factors(.exclude = "q") |>
            add_factors(.exclude = NULL) |>
            ggplot(aes(q, prop, color = outcome)) +
            ggtitle(sprintf("*s* = %s, *u* = %s", season_surv__, u__)) +
            facet_grid( ~ d_yp)
    }
    p <- p +
        geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
        geom_point(aes(shape = outcome)) +
        geom_line(aes(x = as.numeric(q))) +
        scale_y_continuous("Percent of simulations",
                           limits = c(0, 1),
                           breaks = (0:4 * 25) / 100,
                           labels = c("0%", "", "50%", "", "100%")) +
        scale_x_discrete("Between-season determinism (*q*)",
                         labels = c("0", "", "0.5", "", "0.95")) +
        scale_shape_manual(NULL, values = outcome_shapes, drop = FALSE) +
        scale_color_manual(NULL, values = outcome_pal, drop = FALSE) +
        theme(axis.title.x = element_markdown(margin = margin(0,0,0,t=12)),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text.y = element_markdown(),
              plot.title = element_markdown())
    if (season_surv__ != sort(unique(stoch_sim_df$season_surv))[2]) {
        p <- p + theme(axis.title.x = element_blank())
    }
    return(p)
}



# =====================================================*
#           abundances ----
# =====================================================*


abundance_plotter <- function(season_surv__,
                              q__ = NULL,
                              log_trans = TRUE,
                              spp__ = NULL,
                              free_y = FALSE) {

    if (log_trans) {
        .ylimits <- log1p(c(0, 100))
        .ybreaks <- log1p(4^(0:3))
        .ylabels <- c(4^(0:3))
        yvar <- "meanlogY"
        bvar <- "meanlogB"
    } else {
        .ylimits <- c(0, 100)
        .ybreaks <- 0:4 * 25
        .ylabels <- c("0", "", "50", "", "100")
        yvar <- "meanY"
        bvar <- "meanB"
    }

    .facet_fxn <- \(...) facet_grid(scales = "fixed", ...)

    if (free_y) {
        .facet_fxn <- \(...) facet_wrap(scales = "free_y",
                                        ncol = length(unique(stoch_sim_df$d_yp)),
                                        ...)
        .ylimits <- NULL
        .ybreaks <- waiver()
        .ylabels <- waiver()
    }

    if (is.null(q__)) {
        dd <- stoch_sim_df |>
            filter(season_surv == season_surv__) |>
            select(q, u, d_yp, {{ yvar }}, {{ bvar }}) |>
            pivot_longer({{ yvar }}:{{ bvar }}, names_to = "species") |>
            mutate(species = str_sub(species, nchar(species), nchar(species))) |>
            add_factors()
        dds <- dd |>
            group_by(q, u, d_yp, species) |>
            summarize(value = mean(value), .groups = "drop")
        if (!is.null(spp__)) {
            dd <- dd |> filter(species == spp__)
            dds <- dds |> filter(species == spp__)
        }
        p <- dd |>
            ggplot(aes(u, value, color = species)) +
            ggtitle(sprintf("*s* = %s", season_surv__)) +
            .facet_fxn(q ~ d_yp, labeller = labeller(q = q_labeller,
                                                     d_yp = label_value))
    } else {
        dd <- stoch_sim_df |>
            filter(season_surv == season_surv__,
                   q == q__) |>
            select(u, d_yp, {{ yvar }}, {{ bvar }}) |>
            pivot_longer({{ yvar }}:{{ bvar }}, names_to = "species") |>
            mutate(species = str_sub(species, nchar(species), nchar(species))) |>
            add_factors()
        dds <- dd |>
            group_by(u, d_yp, species) |>
            summarize(value = mean(value), .groups = "drop")
        if (!is.null(spp__)) {
            dd <- dd |> filter(species == spp__)
            dds <- dds |> filter(species == spp__)
        }
        p <- dd |>
            ggplot(aes(u, value, color = species)) +
            ggtitle(sprintf("*s* = %s, *q* = %s", season_surv__, q__)) +
            .facet_fxn( ~ d_yp)
    }

    if (!free_y) p <- p +
        geom_hline(yintercept = 0, linetype = 1, color = "gray80")

    p <- p +
        geom_point(aes(color = species), alpha = 0.1, shape = 1) +
        geom_line(data = dds, linewidth = 1) +
        scale_y_continuous("Microbial abundance", limits = .ylimits,
                           breaks = .ybreaks, labels = .ylabels) +
        scale_x_continuous("Strength of microbe--pollinator effect (*u*)",
                           breaks = seq(0, 10, 2.5),
                           labels = c("0", "", "5", "", "10")) +
        scale_shape_manual(NULL, values = c(0, 2), drop = FALSE) +
        scale_linetype_manual(NULL, values = c(1, 1), drop = FALSE) +
        scale_color_manual(NULL, values = spp_pal, drop = FALSE,
                           aesthetics = c("color", "fill")) +
        theme(axis.title.x = element_markdown(margin = margin(0,0,0,t=12)),
              strip.text = element_markdown(),
              plot.title = element_markdown())

    if (season_surv__ != sort(unique(stoch_sim_df$season_surv))[2]) {
        p <- p + theme(axis.title.x = element_blank())
    }

    return(p)
}


abundance_plotter2 <- function(season_surv__,
                               u__ = c(0, 1, 4, 7, 10),
                               log_trans = TRUE,
                               spp__ = NULL,
                               free_y = FALSE) {

    if (log_trans) {
        .ylimits <- log1p(c(0, 100))
        .ybreaks <- log1p(4^(0:3))
        .ylabels <- c(4^(0:3))
        yvar <- "meanlogY"
        bvar <- "meanlogB"
    } else {
        .ylimits <- c(0, 100)
        .ybreaks <- 0:4 * 25
        .ylabels <- c("0", "", "50", "", "100")
        yvar <- "meanY"
        bvar <- "meanB"
    }

    .facet_fxn <- \(...) facet_grid(scales = "fixed", ...)

    if (free_y) {
        .facet_fxn <- \(...) facet_wrap(scales = "free_y",
                                        ncol = length(unique(stoch_sim_df$d_yp)),
                                        ...)
        .ylimits <- NULL
        .ybreaks <- waiver()
        .ylabels <- waiver()
    }

    if (length(u__) > 1) {
        dd <- stoch_sim_df |>
            filter(season_surv == season_surv__) |>
            filter(u %in% u__) |>
            select(q, u, d_yp, {{ yvar }}, {{ bvar }}) |>
            pivot_longer({{ yvar }}:{{ bvar }}, names_to = "species") |>
            mutate(species = str_sub(species, nchar(species), nchar(species))) |>
            add_factors(.exclude = NULL)
        dds <- dd |>
            group_by(q, u, d_yp, species) |>
            summarize(value = mean(value), .groups = "drop")
        if (!is.null(spp__)) {
            dd <- dd |> filter(species == spp__)
            dds <- dds |> filter(species == spp__)
        }
        p <- dd |>
            ggplot(aes(q, value, color = species)) +
            ggtitle(sprintf("*s* = %s", season_surv__)) +
            .facet_fxn(u ~ d_yp, labeller = labeller(u = u_labeller,
                                                     d_yp = label_value))
    } else {
        dd <- stoch_sim_df |>
            filter(season_surv == season_surv__,
                   u == u__) |>
            select(q, d_yp, {{ yvar }}, {{ bvar }}) |>
            pivot_longer({{ yvar }}:{{ bvar }}, names_to = "species") |>
            mutate(species = str_sub(species, nchar(species), nchar(species))) |>
            # add_factors(.exclude = "q")
            add_factors(.exclude = NULL)
        dds <- dd |>
            group_by(q, d_yp, species) |>
            summarize(value = mean(value), .groups = "drop")
        if (!is.null(spp__)) {
            dd <- dd |> filter(species == spp__)
            dds <- dds |> filter(species == spp__)
        }
        p <- dd |>
            ggplot(aes(q, value, color = species)) +
            ggtitle(sprintf("*s* = %s, *u* = %s", season_surv__, u__)) +
            .facet_fxn( ~ d_yp)
    }

    if (!free_y) p <- p +
        geom_hline(yintercept = 0, linetype = 1, color = "gray80")

    p <- p +
        geom_point(aes(color = species), alpha = 0.1, shape = 1) +
        geom_line(data = dds, aes(x = as.numeric(q)), linewidth = 1) +
        scale_y_continuous("Microbial abundance", limits = .ylimits,
                           breaks = .ybreaks, labels = .ylabels) +
        scale_x_discrete("Between-season determinism (*q*)",
                         labels = c("0", "", "0.5", "", "0.95")) +
        scale_shape_manual(NULL, values = c(0, 2), drop = FALSE) +
        scale_linetype_manual(NULL, values = c(1, 1), drop = FALSE) +
        scale_color_manual(NULL, values = spp_pal, drop = FALSE,
                           aesthetics = c("color", "fill")) +
        theme(axis.title.x = element_markdown(margin = margin(0,0,0,t=12)),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_markdown(),
              plot.title = element_markdown())
    if (season_surv__ != sort(unique(stoch_sim_df$season_surv))[2]) {
        p <- p + theme(axis.title.x = element_blank())
    }

    return(p)
}








# =====================================================*
#           gain/loss sims ----
# =====================================================*



if (! file.exists(gain_loss_stoch_sims_file)) {

    # Takes ~12 min with 6 threads
    set.seed(378929239)
    gl_sim_df <- crossing(.u = c(0, 1, 4, 7, 10),
                          .d_yp = d_yp__,
                          .season_surv = sort(unique(stoch_sim_df$season_surv)),
                          .q = sort(unique(stoch_sim_df$q))) |>
        # Don't do this in parallel bc plant_metacomm_stoch is already
        # doing that
        pmap_dfr(\(.u, .d_yp, .season_surv, .q) {
            plant_metacomm_stoch(np = 100L, u = .u, d_yp = .d_yp,
                                 q = .q,
                                 no_immig = TRUE,
                                 n_sigma = 100,
                                 season_surv = .season_surv,
                                 save_every = 150,
                                 max_t = 3000) |>
                select(-P) |>
                pivot_longer(Y:B, names_to = "species", values_to = "density") |>
                mutate(season = as.integer((t - 0.1) %/% 150 + 1)) |>
                mutate(occup = density > 0) |>
                group_by(rep, p, species) |>
                # this season is occupied and next one is unoccupied:
                mutate(loss = occup & lead(!occup)) |>
                # this season is unoccupied and next one is occupied:
                mutate(gain = !occup & lead(occup)) |>
                ungroup() |>
                filter(season < max(season)) |>
                group_by(rep, season, species) |>
                # total landscape-wide loss / gain and abundance
                summarize(across(c(density, occup, loss, gain), sum),
                          .groups = "drop") |>
                mutate(u = .u, d_yp = .d_yp, season_surv = .season_surv,
                       q = .q) |>
                add_factors(.exclude = NULL) |>
                select(u, d_yp, season_surv, q, everything())
        }, .progress = TRUE)

    write_rds(gl_sim_df, gain_loss_stoch_sims_file)

} else {

    gl_sim_df <- read_rds(gain_loss_stoch_sims_file)

}


#'
#' Average proportion gained or lost (or net = gained - lost)
#'
prop_gl_sims <- map(c("gain", "loss", "net"), \(type__) {

    .y_axis <- scale_y_continuous(sprintf("Average proportion %s",
                                          case_when(type__ == "gain" ~ "gained",
                                                    type__ == "loss" ~ "lost",
                                                    TRUE ~ type__)),
                                  limits = c(0, 1.1), breaks = 0:4/4,
                                  labels = c("0", "", "0.5", "", "1.0"))
    if (type__ == "net") {
        .y_axis$limits <- NULL
        .y_axis$breaks <- waiver()
        .y_axis$labels <- waiver()
    }
    gl_sim_df |>
        filter(season_surv == median(season_surv)) |>
        filter(occup > 0) |>
        mutate(net = gain - loss) |>
        mutate(across(c(net, gain, loss), \(x) x / occup)) |>
        # mean across reps and seasons:
        group_by(u, d_yp, q, species) |>
        summarize(across(c(net, gain, loss), mean), .groups = "drop") |>
        mutate(q = as.numeric(paste(q))) |>
        ggplot(aes(q, .data[[type__]], color = species)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_point(size = 2) +
        geom_line(linewidth = 1) +
        scale_x_continuous("Between-season determinism (*q*)",
                           breaks = as.numeric(levels(gl_sim_df$q)),
                           labels = c("0", "", "0.5", "", "0.95")) +
        .y_axis +
        facet_grid(u ~ d_yp,
                   labeller = labeller(u = u_labeller, d_yp = label_value)) +
        scale_color_manual(NULL, values = spp_pal) +
        theme(axis.title.x = element_markdown(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_markdown())
})

do.call(wrap_plots, prop_gl_sims[1:3]) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold", vjust = -10))



# REALLY LEFT OFF ----
# gl_plot <-

map(1:2, \(i) {

})

dd <- gl_sim_df |>
    filter(season_surv == median(season_surv)) |>
    filter(occup > 0) |>
    mutate(net = gain - loss) |>
    mutate(across(c(net, gain, loss), \(x) x / occup)) |>
    # mean across reps and seasons:
    group_by(u, d_yp, q, species) |>
    summarize(across(c(net, gain, loss), mean), .groups = "drop") |>
    mutate(q = as.numeric(paste(q)))

if (dd == 1) {
    dd <- dd |>
        select(-net) |>
        pivot_longer(gain:loss)
    .y_axis <- scale_y_continuous("Average proportion gained or lost",
                                  limits = c(0, 1.1), breaks = 0:4/4,
                                  labels = c("0", "", "0.5", "", "1.0"))
} else {
    dd <- dd |>
        mutate(name = "net", value = net)
    .y_axis <- scale_y_continuous("Average proportion net-gained")
}

dd |>
    ggplot(aes(q, value, color = species)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
    geom_point(aes(shape = name), size = 1) +
    geom_line(aes(linetype = name), linewidth = 1) +
    scale_x_continuous("Between-season determinism (*q*)",
                       breaks = as.numeric(levels(gl_sim_df$q)),
                       labels = c("0", "", "0.5", "", "0.95")) +
    .y_axis +
    facet_grid(u ~ d_yp,
               labeller = labeller(u = u_labeller, d_yp = label_value)) +
    scale_color_manual(NULL, values = spp_pal) +
    scale_shape_manual(values = c(0, 2)) +
    scale_linetype_manual(values = c(net = "solid", gain = "solid",
                                     loss = "22")) +
    theme(axis.title.x = element_markdown(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.title = element_blank(),
          strip.text = element_markdown())






# gl_sim_df |>
#     filter(season_surv == median(season_surv)) |>
#     filter(occup > 0) |>
#     mutate(net = (gain - loss) / occup) |>
#     select(u, d_yp, q, rep, season, species, net) |>
#     pivot_wider(names_from = species, values_from = net) |>
#     filter(!is.na(yeast), !is.na(bacteria)) |>
#     mutate(diff = yeast - bacteria,
#            diff2 = ifelse(d_yp == "bacteria only",
#                           yeast - bacteria, bacteria - yeast)) |>
#     # mean across reps and seasons:
#     group_by(u, d_yp, q) |>
#     summarize(across(starts_with("diff"), mean), .groups = "drop") |>
#     mutate(q = as.numeric(paste(q))) |>
#     ggplot(aes(q, diff2)) +
#     geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
#     geom_point(size = 2) +
#     geom_line(linewidth = 1) +
#     scale_x_continuous("Between-season determinism (*q*)",
#                        breaks = as.numeric(levels(gl_sim_df$q)),
#                        labels = c("0", "", "0.5", "", "0.95")) +
#     ylab("rare - common") +
#     # ylab("yeast - bacteria") +
#     # facet_grid(u ~ d_yp,
#     facet_wrap(u ~ d_yp,
#                scales = "free_y", ncol = length(unique(stoch_sim_df$d_yp)),
#                labeller = labeller(u = u_labeller, d_yp = label_value)) +
#     scale_color_manual(values = spp_pal, guide = "none") +
#     theme(axis.title.x = element_markdown(),
#           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           strip.text = element_markdown(),
#           plot.title = element_markdown())







# |>
#     group_by(u, d_yp, q, rep, species) |>
#     summarize(sumnet = sum(net), .groups = "drop") |>
#     pivot_wider(names_from = species, values_from = sumnet) |>
#     mutate(diff = yeast - bacteria,
#            q = as.numeric(paste(q)))
# dds <- dd |>
#     group_by(u, d_yp, q) |>
#     summarize(diff = mean(diff), .groups = "drop")
#
# dd |>
#     ggplot(aes(q, diff)) +
#     ggtitle(sprintf("*s* = %s", season_surv__)) +
#     geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
#     geom_point(alpha = 0.1, shape = 1) +
#     geom_line(data = dds, linewidth = 1) +
#     scale_y_continuous("Net gain yeast - net gain bacteria",
#                        limits = .ylimits,
#                        breaks = -1:2 * 50) +
#     scale_x_continuous("Between-season determinism (*q*)",
#                        breaks = as.numeric(levels(gl_sim_df$q)),
#                        labels = c("0", "", "0.5", "", "0.95")) +
#     facet_grid(u ~ d_yp,
#                labeller = labeller(u = u_labeller, d_yp = label_value)) +
#     scale_color_manual(values = spp_pal, guide = "none") +
#     theme(axis.title.x = element_markdown(),
#           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           strip.text.y = element_markdown(),
#           plot.title = element_markdown())






#' Plot difference in net gain/loss for yeast - net gain/loss for bacteria
#' `ylims_by_ss` is whether to have separate y limits by `season_surv__`
diff_net_plotter <- function(season_surv__,
                             ylims_by_ss = TRUE) {

    .ylimits <- NULL

    if (!ylims_by_ss) {
        .ylimits <- gl_sim_df |>
            mutate(net = gain - loss) |>
            group_by(u, d_yp, season_surv, q, rep, species) |>
            summarize(sumnet = sum(net), .groups = "drop") |>
            pivot_wider(names_from = species, values_from = sumnet) |>
            mutate(diff = yeast - bacteria) |>
            getElement("diff") |>
            (\(x) range(x) + c(-1L, 1L))()
    }


    dd <- gl_sim_df |>
        filter(season_surv == season_surv__) |>
        mutate(net = gain - loss) |>
        group_by(u, d_yp, q, rep, species) |>
        summarize(sumnet = sum(net), .groups = "drop") |>
        pivot_wider(names_from = species, values_from = sumnet) |>
        mutate(diff = yeast - bacteria,
               q = as.numeric(paste(q)))
    dds <- dd |>
        group_by(u, d_yp, q) |>
        summarize(diff = mean(diff), .groups = "drop")

    dd |>
        ggplot(aes(q, diff)) +
        ggtitle(sprintf("*s* = %s", season_surv__)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_point(alpha = 0.1, shape = 1) +
        geom_line(data = dds, linewidth = 1) +
        scale_y_continuous("Net gain yeast - net gain bacteria",
                           limits = .ylimits,
                           breaks = -1:2 * 50) +
        scale_x_continuous("Between-season determinism (*q*)",
                           breaks = as.numeric(levels(gl_sim_df$q)),
                           labels = c("0", "", "0.5", "", "0.95")) +
        facet_grid(u ~ d_yp,
                   labeller = labeller(u = u_labeller, d_yp = label_value)) +
        scale_color_manual(values = spp_pal, guide = "none") +
        theme(axis.title.x = element_markdown(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text.y = element_markdown(),
              plot.title = element_markdown())
}


#' Plot net gain/loss for yeast and bacteria separately
net_plotter <- function(season_surv__,
                        ylims_by_ss = TRUE) {

    .ylimits <- NULL

    if (!ylims_by_ss) {
        .ylimits <- gl_sim_df |>
            mutate(net = gain - loss) |>
            group_by(u, d_yp, season_surv, q, rep, species) |>
            summarize(sumnet = sum(net), .groups = "drop") |>
            getElement("sumnet") |>
            (\(x) range(x) + c(-1L, 1L))()
    }


    dd <- gl_sim_df |>
        filter(season_surv == season_surv__) |>
        mutate(net = gain - loss) |>
        group_by(u, d_yp, q, rep, species) |>
        summarize(sumnet = sum(net), .groups = "drop") |>
        mutate(q = as.numeric(paste(q)))
    dds <- dd |>
        group_by(u, d_yp, q, species) |>
        summarize(sumnet = mean(sumnet), .groups = "drop")

    dd |>
        ggplot(aes(q, sumnet, color = species)) +
        ggtitle(sprintf("*s* = %s", season_surv__)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_point(alpha = 0.1, shape = 1) +
        geom_line(data = dds, linewidth = 1) +
        scale_y_continuous("Summed net occupancy gains",
                           limits = .ylimits,
                           breaks = -1:2 * 50) +
        scale_x_continuous("Between-season determinism (*q*)",
                           breaks = as.numeric(levels(gl_sim_df$q)),
                           labels = c("0", "", "0.5", "", "0.95")) +
        facet_grid(u ~ d_yp,
                   labeller = labeller(u = u_labeller, d_yp = label_value)) +
        scale_color_manual(values = spp_pal, guide = "none") +
        theme(axis.title.x = element_markdown(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text.y = element_markdown(),
              plot.title = element_markdown())
}





gl_plotter <- function(season_surv__,
                       type__,
                       ylims_by_ss = TRUE) {

    type__ <- match.arg(type__, c("gain", "loss"))

    if (!ylims_by_ss) {
        .ylimits <- gl_sim_df |>
            group_by(u, d_yp, season_surv, q, rep, species) |>
            summarize(loss = sum(loss), gain = sum(gain), .groups = "drop") |>
            pivot_longer(loss:gain) |>
            group_by(u, d_yp, season_surv, q, species, name) |>
            summarize(value = mean(value), .groups = "drop") |>
            getElement("value") |>
            (\(x) c(0L, max(x) + 1L))()
    } else {
        .ylimits <- gl_sim_df |>
            filter(season_surv == season_surv__) |>
            group_by(u, d_yp, q, rep, species) |>
            summarize(loss = sum(loss), gain = sum(gain), .groups = "drop") |>
            pivot_longer(loss:gain) |>
            group_by(u, d_yp, q, species, name) |>
            summarize(value = mean(value), .groups = "drop") |>
            getElement("value") |>
            (\(x) c(0L, max(x) + 1L))()
    }

    gl_sim_df |>
        filter(season_surv == season_surv__) |>
        group_by(u, d_yp, q, rep, species) |>
        summarize({{ type__ }} := sum(.data[[type__]]), .groups = "drop") |>
        group_by(u, d_yp, q, species) |>
        summarize(value = mean(.data[[type__]]), .groups = "drop") |>
        mutate(q = as.numeric(paste(q))) |>
        ggplot(aes(q, value, color = species)) +
        ggtitle(sprintf("*s* = %s", season_surv__)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_point() +
        geom_line(linewidth = 1) +
        scale_y_continuous(sprintf("Total %s over all seasons",
                                   ifelse(type__ == "gain", "gains", "losses")),
                           limits = .ylimits, breaks = 0:2 * 250) +
        scale_x_continuous("Between-season determinism (*q*)",
                           breaks = as.numeric(levels(gl_sim_df$q)),
                           labels = c("0", "", "0.5", "", "0.95")) +
        facet_grid(u ~ d_yp,
                   labeller = labeller(u = u_labeller, d_yp = label_value)) +
        scale_color_manual(values = spp_pal, guide = "none") +
        theme(axis.title.x = element_markdown(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_markdown(),
              plot.title = element_markdown())


    # gl_sim_df |>
    #     filter(season_surv == season_surv__) |>
    #     group_by(u, d_yp, q, rep, species) |>
    #     summarize(loss = sum(loss), gain = sum(gain), .groups = "drop") |>
    #     pivot_longer(loss:gain) |>
    #     mutate(q = as.numeric(paste(q))) |>
    #     group_by(u, d_yp, q, species, name) |>
    #     summarize(value = mean(value), .groups = "drop") |>
    #     ggplot(aes(q, value, color = species)) +
    #     ggtitle(sprintf("*s* = %s", season_surv__)) +
    #     geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
    #     geom_point(aes(shape = name)) +
    #     geom_line(aes(linetype = name), linewidth = 1) +
    #     scale_y_continuous("Total gains or losses over 10 seasons",
    #                        limits = .ylimits,
    #                        breaks = 0:2 * 250) +
    #     scale_x_continuous("Between-season determinism (*q*)",
    #                        breaks = as.numeric(levels(gl_sim_df$q)),
    #                        labels = c("0", "", "0.5", "", "0.95")) +
    #     facet_grid(u ~ d_yp,
    #                labeller = labeller(u = u_labeller, d_yp = label_value)) +
    #     scale_color_manual(values = spp_pal, guide = "none") +
    #     scale_shape_manual(NULL, values = c(0, 2)) +
    #     scale_linetype_manual(NULL, values = c("solid", "22")) +
    #     theme(axis.title.x = element_markdown(),
    #           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #           strip.text = element_markdown(),
    #           plot.title = element_markdown())
}


(net_plotter(median(gl_sim_df$season_surv)) +
        diff_net_plotter(median(gl_sim_df$season_surv)) &
        theme(plot.title = element_blank())) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))


# map(sort(unique(gl_sim_df$season_surv)), gl_plotter) |>
#     do.call(what = wrap_plots) +
#     plot_layout(guides = "collect", nrow = 1)

gl_plotter(median(gl_sim_df$season_surv))


# LEFT OFF ----

#' What to include in figure (all are only for s = 0.02):
#'   1. u vs outcomes (just for `q = 0.5`)
#'   2. u vs microbial abundance (just for `q = 0.5`)
#'   3. Something showing how non-random between seasons makes coexistence
#'      less likely at low s and more likely at high s?
#'


outcome_plotter(median(unique(stoch_sim_df$season_surv)), q__ = 0.5) /
    abundance_plotter(median(unique(stoch_sim_df$season_surv)), q__ = 0.5) &
    theme(plot.title = element_blank())


# For supplement

supp_q_out_abund_p <-
    (outcome_plotter(median(unique(stoch_sim_df$season_surv))) +
         theme(plot.title = element_blank())) +
    (abundance_plotter(median(unique(stoch_sim_df$season_surv))) +
         theme(plot.title = element_blank())) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))


supp_u_out_abund_p <-
    (outcome_plotter2(median(unique(stoch_sim_df$season_surv))) +
    theme(plot.title = element_blank())) +
    (abundance_plotter2(median(unique(stoch_sim_df$season_surv))) +
         theme(plot.title = element_blank())) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))

# save_plot("_figures/supp-q-out-abund.pdf", supp_q_out_abund_p, 9, 6)
# save_plot("_figures/supp-u-out-abund.pdf", supp_u_out_abund_p, 9, 6)



