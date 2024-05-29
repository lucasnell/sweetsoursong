

library(sweetsoursong)
library(tidyverse)
library(parallel)

options("mc.cores" = max(1L, detectCores()-2))

outcome_pal <- c("coexist" = "#993399", # "#4477AA",
                 "coexist - within patch" = "#CC6666",  # "#66CCEE",
                 "yeast only" = "#FFCC33", # "#CCBB44",
                 "bacteria only" = "#333399") # "#AA3377")

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "gray60")



even_run <- function(other_args, np, no_error = FALSE) {

    if (inherits(other_args, "data.frame")) {
        stopifnot(nrow(other_args) == 1L)
        other_args <- as.list(other_args)
    }

    run_fun <- landscape_constantF_stoch_ode
    if (no_error) run_fun <- landscape_constantF_ode

    stopifnot(length(np) == 1 && is.numeric(np) && np %% 1 == 0 && np >= 2)
    np <- as.integer(np)

    args <- list(m = 0.1,
                 # ----------*
                 # within-plant dispersal rates:
                 d_yp = 1.1,  # 1.5,  # pollinator-dependent for yeast (**)
                 d_b0 = 0.3,  # pollinator-independent for bacteria (**)
                 d_bp = 0.4,  # pollinator-dependent for bacteria (**)
                 # -------------*
                 # rates of immigration from non-focal-plant sources:
                 g_yp = 0, # 0.005, # pollinator-dependent for yeast (**)
                 g_b0 = 0, # 0.02, # pollinator-independent for bacteria (**)
                 g_bp = 0, # 0.002,  # pollinator-dependent for bacteria
                 # -------------*
                 # pollinators:
                 L_0 = 1 / np,  # half saturation ratio for P -> dispersal (--)
                 u = 0,      # power for B -> P (??)
                 X = 0,      # attraction for non-focal plants (??)
                 #
                 # (**)  = value from Song et al. (submitted)
                 # (--)  = value different from Song et al. (submitted)
                 # (??) = value should be varied bc I have no clue what to use
                 #
                 sigma_y = 0,
                 sigma_b = 0,
                 dt = 0.1,
                 max_t = 120)

    for (n in names(args)[!names(args) %in% c("max_t","dt","X","u","sigma_y","sigma_b")]) {
        args[[n]] <- rep(args[[n]], np)
    }

    stopifnot(inherits(other_args, "list"))
    stopifnot(! is.null(names(other_args)) && ! any(names(other_args) == ""))

    for (n in names(other_args)) {
        if (! n %in% c("P_total", "u", "X", "max_t","sigma_y","sigma_b") &&
            length(other_args[[n]]) == 1) {
            other_args[[n]] <- rep(other_args[[n]], np)
        }
        args[[n]] <- other_args[[n]]
    }


    args[["Y0"]] <- seq(0.1, 0.4, length.out = np)
    args[["B0"]] <- 0.5 - args[["Y0"]]

    if (no_error) {
        args[["sigma_y"]] <- NULL
        args[["sigma_b"]] <- NULL
    }

    run_df <- do.call(run_fun, args) |>
        as_tibble() |>
        mutate(p = factor(p, levels = 0:(np-1L),
                          labels = paste("patch", 1:np)))

    return(run_df)

}


one_even_plot <- function(x, ymax = NULL, no_labs = FALSE) {
    p <- x |>
        # filter(t %% 1 == 0) |>
        pivot_longer(Y:P, names_to = "type", values_to = "density") |>
        mutate(type = factor(type, levels = c("Y", "B", "P"),
                             labels = c("yeast", "bacteria", "pollinators"))) |>
        ggplot(aes(t, density)) +
        geom_hline(yintercept = 0, linewidth = 1, color = "gray80") +
        geom_line(aes(color = type, linetype = type), linewidth = 1) +
        facet_grid( ~ p) +
        xlab("Time (days)") +
        scale_color_manual(NULL, values = spp_pal) +
        scale_linetype_manual(NULL, values = c("solid", "solid", "24"))
    if (no_labs) p <- p +
            theme(strip.text = element_blank(), strip.background = element_blank(),
                  legend.position = "none", axis.title = element_blank())
    if (!is.null(ymax)) p <- p + ylim(0, ymax)
    return(p)
}


list(u = 0, d_yp = 0.95, sigma_y = 0.1, sigma_b = 0.1) |>
    even_run(np = 5) |>
    one_even_plot()

list(u = 0.25, d_yp = 0.95, max_t = 2000) |>
    even_run(np = 5) |>
    one_even_plot()

