
#'
#' This script creates the following files inside `_figures`:
#'
#'   - `u-stochastic-main.pdf`
#'   - `supp-u-out-abund.pdf`
#'
#' File `u-stochastic-main.pdf` is used inside figure 4 to show outcomes
#' and microbial abundances across values of `u` and `d_yp` for stochastic
#' simulations of 100-plant landscapes when both yeast and bacteria
#' start with the same abundance.
#' This is to show how `u` affects long-term co-occurrence of competitors.
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

if (! file.exists(stoch_rds_files$main)) {

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

    write_rds(stoch_sim_df, stoch_rds_files$main)

} else {

    stoch_sim_df <- read_rds(stoch_rds_files$main)

}





# ============================================================================*
# ============================================================================*
# Functions ----
# ============================================================================*
# ============================================================================*

#'
#' Add the proportion of different outcomes to a data frame.
#'
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

#'
#' Create a plot of outcomes
#'
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

#'
#' Create a plot of abundances.
#'
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







# ============================================================================*
# ============================================================================*
# Create plots ----
# ============================================================================*
# ============================================================================*


#'
#' Main text figure panels with `q = 0.5`:
#'
main_plots <- list(outcome_plotter("u", facet_q = 0.5, facet_u = 4),
                   abundance_plotter("u", facet_q = 0.5, facet_u = 4)) |>
    map(\(p) {
        p + subpanel_theme + theme(panel.spacing.x = unit(1, "lines"))
    }) |>
    do.call(what = wrap_plots) +
    plot_layout(ncol = 1)

if (.write_plots) {
    save_plot("_figures/u-stochastic-main.pdf", main_plots, 4, 3.4)
} else {
    main_plots
}



#'
#' Supplemental figure of how `u` affects outcomes and abundances from
#' evenly distributed starting abundances, across values of  `q`:
#'
supp_u_out_abund_p <-
    (outcome_plotter("u", drop_ext = FALSE) +
         theme(plot.title = element_blank())) +
    (abundance_plotter("u") + theme(plot.title = element_blank())) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))

if (.write_plots) {
    save_plot("_figures/supp-u-out-abund.pdf", supp_u_out_abund_p, 9, 6)
} else {
    supp_u_out_abund_p
}
