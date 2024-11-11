
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
# plot functions ----
# ===========================================================================*
# ===========================================================================*

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
                            add_title = FALSE) {

    x_var <- match.arg(x_var, c("q", "u"))

    # makes objects: x_scale, facet, f_lvls, f_var, f_lvls
    make_plot_objects(x_var, facet_q, facet_u,
                      env = sys.frame(sys.nframe()))


    p <- stoch_sim_df |>
        filter(season_surv == median(season_surv),
               .data[[f_var]] %in% f_lvls) |>
        select(q, u, d_yp, maxY, maxB) |>
        add_outcome_props() |>
        add_factors(.exclude = x_var) |>
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
        summarize(value = mean(value), .groups = "drop")
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





# =====================================================*
#           creating plots ----
# =====================================================*



# LEFT OFF ----

#' What to include in figure (all are only for s = 0.02):
#'   1. u vs outcomes (just for `q = 0.5`)
#'   2. u vs microbial abundance (just for `q = 0.5`)
#'   3. Something showing how non-random between seasons makes coexistence
#'      less likely at low s and more likely at high s?
#'


outcome_plotter("u", facet_q = 0.5) /
    abundance_plotter("u", facet_q = 0.5) &
    theme(plot.title = element_blank())

outcome_plotter("q", facet_u = 4) /
    abundance_plotter("q", facet_u = 4) &
    theme(plot.title = element_blank())


# For supplement

supp_u_out_abund_p <-
    (outcome_plotter("u") + theme(plot.title = element_blank())) +
    (abundance_plotter("u") + theme(plot.title = element_blank())) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))


supp_q_out_abund_p <-
    (outcome_plotter("q") + theme(plot.title = element_blank())) +
    (abundance_plotter("q") + theme(plot.title = element_blank())) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))

# save_plot("_figures/supp-u-out-abund.pdf", supp_u_out_abund_p, 9, 6)
# save_plot("_figures/supp-q-out-abund.pdf", supp_q_out_abund_p, 9, 6)




#'
#' Average proportion (1) gained or lost and (2) net (gained - lost)
#'
gl_plots <- map(1:2, \(i) {
    dd <- gl_sim_df |>
        filter(season_surv == median(season_surv)) |>
        filter(occup > 0) |>
        mutate(net = gain - loss) |>
        mutate(across(c(net, gain, loss), \(x) x / occup)) |>
        # mean across reps and seasons:
        group_by(u, d_yp, q, species) |>
        summarize(across(c(net, gain, loss), mean), .groups = "drop") |>
        mutate(q = as.numeric(paste(q)))
    if (i == 1) {
        dd <- dd |>
            select(-net) |>
            pivot_longer(gain:loss) |>
            mutate(name = factor(name, levels = c("gain", "loss", "net")))
        .y_axis <- scale_y_continuous("Average proportion gained or lost",
                                      limits = c(0, 1.1), breaks = 0:4/4,
                                      labels = c("0", "", "0.5", "", "1.0"))
        .guide <- guide_legend()
    } else {
        dd <- dd |>
            mutate(name = factor("net", levels = c("gain", "loss", "net")),
                   value = net)
        .y_axis <- scale_y_continuous("Average proportion net-gained")
        .guide <- "none"
    }
    u_labeller <- \(x) lapply(x, \(z) paste('*u* = ', z))
    dd |>
        ggplot(aes(q, value, color = species)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_point(aes(shape = name), size = 2, show.legend = TRUE) +
        geom_line(aes(linetype = name), linewidth = 1, show.legend = TRUE) +
        scale_x_continuous("Between-season determinism (*q*)",
                           breaks = as.numeric(levels(gl_sim_df$q)),
                           labels = c("0", "", "0.5", "", "0.95")) +
        .y_axis +
        facet_grid(u ~ d_yp,
                   labeller = labeller(u = u_labeller,
                                       d_yp = label_value)) +
        scale_color_manual(NULL, values = spp_pal) +
        scale_shape_manual(values = c(gain = 0, loss = 2, net = 19),
                           drop = FALSE) +
        scale_linetype_manual(values = c(gain = "solid", loss = "22",
                                         net = "solid"),
                              drop = FALSE) +
        theme(axis.title.x = element_markdown(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.title = element_blank(),
              strip.text = element_markdown())
})

supp_q_gain_loss_p <- do.call(wrap_plots, gl_plots) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))

# save_plot("_figures/supp-q-gain-loss.pdf", supp_q_gain_loss_p, 9, 6)

