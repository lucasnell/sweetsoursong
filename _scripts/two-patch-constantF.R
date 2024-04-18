

library(sweetsoursong)
library(tidyverse)
library(parallel)

options("mc.cores" = max(1L, detectCores()-2))

outcome_pal <- c("coexist" = "#CB3694", # "#4477AA",
                 "coexist - within patch" = "#CC6666",  # "#66CCEE",
                 "yeast only" = "#FFCC33", # "#CCBB44",
                 "bacteria only" = "#333399") # "#AA3377")

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "gray60")



high_low_run <- function(other_args) {

    if (inherits(other_args, "data.frame")) {
        stopifnot(nrow(other_args) == 1L)
        other_args <- as.list(other_args)
    }

    args <- list(dt = 0.1,
                 max_t = 100,
                 m = rep(0.1, 2))

    stopifnot(inherits(other_args, "list"))
    stopifnot(! is.null(names(other_args)) && ! any(names(other_args) == ""))

    if ("max_t" %in% names(other_args)) {
        args[["max_t"]] <- other_args[["max_t"]]
        other_args[["max_t"]] <- NULL
    }

    stopifnot(! any(names(other_args) %in% names(args)))

    for (n in names(other_args)) {
        if (n %in% c("P_total", "S_0", "X")) next
        if (length(other_args[[n]]) == 1) other_args[[n]] <- rep(other_args[[n]], 2)
    }

    # Below, I made the low-density one a bit different to avoid weirdness
    # with exactly equal starting values
    all_runs_args <-list("high B" = list(Y0 = c(0.11, 0.1),  B0 = c(0.5, 0.5)),
                         "high Y" = list(Y0 = c(0.5, 0.5), B0 = c(0.1, 0.11)),
                         "one of each" = list(Y0 = c(0.5, 0.1), B0 = c(0.1, 0.5)))

    all_runs <- all_runs_args |>
        imap_dfr(\(x, n) {
            c(args, other_args, x) |>
                do.call(what = landscape_constantF_ode) |>
                as_tibble() |>
                mutate(run = n)
        }) |>
        select(run, everything()) |>
        mutate(run = factor(run, levels = names(all_runs_args)),
               p = factor(p, levels = 0:1, labels = c("patch 1", "patch 2")))

    return(all_runs)

}


tibble(# ----------
       # within-plant dispersal rates:
       d_yp = 1.0, # 1.5,  # pollinator-dependent for yeast (**)
       d_b0 = 0.3,  # pollinator-independent for bacteria (**)
       d_bp = 0.4,  # pollinator-dependent for bacteria (**)
       # -------------
       # rates of immigration from non-focal-plant sources:
       g_yp = 0, # 0.005, # pollinator-dependent for yeast (**)
       g_b0 = 0, # 0.02, # pollinator-independent for bacteria (**)
       g_bp = 0, # 0.002,  # pollinator-dependent for bacteria
       # -------------
       # pollinators:
       L_0 = 0.5,          # half saturation ratio for P -> dispersal (**)
       S_0 = 10,              # strength of B -> P (??)
       X = 0,              # attraction for nearby flowers other than focal plant (??)

       #
       # (**)  = value from Song et al. (submitted)
       # (??) = value should be varied bc I have no clue what to use
       #
       max_t = 500
) |>
    high_low_run() |>
    (\(x) {
        print(group_by(x, run) |>
                  filter(t == max(t)) |>
                  summarize(Y = mean(Y),
                            B = mean(B),
                            N = 1 - Y - B,
                            .groups = "drop"))
        return(x)
    })() |>
    filter(t %% 1 == 0) |>
    pivot_longer(Y:P, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B", "P"),
                         labels = c("yeast", "bacteria", "pollinators"))) |>
    ggplot(aes(t, density)) +
    geom_hline(yintercept = 0, linewidth = 1, color = "gray80") +
    geom_line(aes(color = type, linetype = type), linewidth = 1) +
    facet_grid(run ~ p) +
    xlab("Time (days)") +
    scale_color_manual(NULL, values = spp_pal) +
    scale_linetype_manual(NULL, values = c("solid", "solid", "24")) +
    # theme(strip.text = element_blank(), strip.background = element_blank(),
    #       legend.position = "none") +
    NULL





equil_run <- function(other_args, .n_plants = 2L, .precision = 1e-4) {

    stopifnot(inherits(other_args, "list"))
    stopifnot(all(map_int(other_args, length) == 1))

    args <- list(dt = 0.1,
                 max_t = 1000,
                 m = rep(0.1, .n_plants),
                 # -------------
                 # rates of immigration from non-focal-plant sources:
                 g_yp = rep(0,.n_plants), # 0.005, # pollinator-dependent for yeast (**)
                 g_b0 = rep(0,.n_plants), # 0.02, # pollinator-independent for bacteria (**)
                 g_bp = rep(0,.n_plants), # 0.002,  # pollinator-dependent for bacteria
                 # -------------
                 # pollinators:
                 X = 0               # attraction for nearby flowers other than focal plant (??)
    )

    stopifnot(inherits(other_args, "list"))
    stopifnot(! is.null(names(other_args)) && ! any(names(other_args) == ""))

    if ("max_t" %in% names(other_args)) {
        args[["max_t"]] <- other_args[["max_t"]]
        other_args[["max_t"]] <- NULL
    }

    stopifnot(! any(names(other_args) %in% names(args)))

    for (n in names(other_args)) {
        if (n %in% c("P_total", "S_0", "X")) next
        if (length(other_args[[n]]) == 1) other_args[[n]] <- rep(other_args[[n]], .n_plants)
    }

    # Below, I made the low-density one a bit different to avoid weirdness
    # with exactly equal starting values
    all_runs_args <-list("high B" = list(Y0 = c(0.11, 0.1),  B0 = c(0.5, 0.5)),
                         "high Y" = list(Y0 = c(0.5, 0.5), B0 = c(0.1, 0.11)),
                         "one of each" = list(Y0 = c(0.5, 0.1), B0 = c(0.1, 0.5)))


    all_runs <- all_runs_args |>
        imap_dfr(\(x, n) {
            o <- c(args, other_args, x) |>
                do.call(what = landscape_constantF_ode) |>
                as_tibble() |>
                filter(t >= max(t) - 0.1)
            for (spp in c("Y","B")) {
                Nt1 <- o[[spp]][o$t == max(o$t)]
                Nt0 <- o[[spp]][o$t == min(o$t)]
                dyb <- max(abs(Nt1 - Nt0))
                if (dyb > .precision) stop("More time points needed for ", n)
            }
            o |>
                filter(t == max(t)) |>
                mutate(p = p + 1L, run = n) |>
                select(-t) |>
                pivot_wider(names_from = p, values_from = c(Y, B, P),
                            names_sep = "")
        }) |>
        select(run, everything()) |>
        mutate(run = factor(run, levels = names(all_runs_args)))

    for (n in names(other_args)) {
        all_runs[[n]] <- other_args[[n]][[1]]
    }

    return(all_runs)

}





if (! file.exists("_data/two-patch-constF-equil.rds")) {

    #' 2700 with 6 threads takes ~2 min
    equil_df <- crossing(# ----------
                         # within-plant dispersal rates:
                         d_yp = 1.1 + -4:4/10, # 1.5,  # pollinator-dependent for yeast (**)
                         d_b0 = 0.3 + -2:2/10,  # pollinator-independent for bacteria (**)
                         d_bp = 0.4 + -3:1/10,  # pollinator-dependent for bacteria (**)
                         # -------------
                         # pollinators:
                         L_0 = c(0.1, 0.5, 0.9),          # for P -> dispersal (**)
                         S_0 = c(0, 1, 5, 10),              # strength of B -> P (??)
                         max_t = 10e3L) |>
        (\(x) split(x, 1:nrow(x)))() |>
        map(as.list) |>
        mclapply(equil_run) |>
        do.call(what = bind_rows)

    write_rds(equil_df, "_data/two-patch-constF-equil.rds")

} else {

    equil_df <- read_rds("_data/two-patch-constF-equil.rds")

}


outcome_plot_prep <- function(x) {
    x |>
        mutate(Yprop = (Y1 + Y2) / (Y1 + Y2 + B1 + B2),
               S_0 = factor(S_0, levels = sort(unique(S_0)),
                          labels = sprintf("S_0 = %i", sort(unique(S_0)))),
               outcome = case_when(Yprop < 1e-3 ~ "bacteria only",
                                   Yprop > (1-1e-3) ~ "yeast only",
                                   Y1 >= 1e-3 & Y1 <= (1-1e-3) &
                                       Y2 >= 1e-3 & Y2 <= (1-1e-3) ~ "coexist - within patch",
                                   TRUE ~ "coexist"))
}



equil_df |>
    filter(L_0 == 0.5, d_bp == 0.4) |>
    outcome_plot_prep() |>
    ggplot(aes(d_yp, d_b0, fill = outcome)) +
    geom_raster() +
    scale_fill_manual(values = outcome_pal) +
    xlab("Yeast pollinator-dependent dispersal") +
    ylab("Bacteria pollinator-independent dispersal") +
    facet_grid(run ~ S_0)


equil_df |>
    filter(L_0 == 0.5, d_b0 == 0.3) |>
    outcome_plot_prep() |>
    ggplot(aes(d_yp, d_bp, fill = outcome)) +
    geom_raster() +
    scale_fill_manual(values = outcome_pal) +
    xlab("Yeast pollinator-dependent dispersal") +
    ylab("Bacteria pollinator-dependent dispersal") +
    facet_grid(run ~ S_0)



# For Yosemite talk
equil_df |>
    filter(L_0 == 0.5, d_bp == 0.4) |>
    filter(run == "one of each", S_0 %in% c(0, 1, 10)) |>
    outcome_plot_prep() |>
    ggplot(aes(d_yp, d_b0, fill = outcome)) +
    geom_raster() +
    scale_fill_manual(values = outcome_pal) +
    xlab("Yeast pollinator-dependent dispersal") +
    ylab("Bacteria pollinator-independent dispersal") +
    facet_wrap(~ S_0, nrow = 1) +
    theme(legend.position = "none", axis.title = element_blank(),
          strip.text = element_blank(), strip.background = element_blank())






