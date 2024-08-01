

library(sweetsoursong)
library(tidyverse)
library(patchwork)
library(RcppParallel)
setThreadOptions(numThreads = max(defaultNumThreads() - 2L, 1L))

if (file.exists(".Rprofile")) source(".Rprofile")

stoch_sims_file <- "_data/constantF-stoch-sims.rds"
open_stoch_sims_file <- "_data/constantF-stoch-open-sims.rds"


outcome_pal <- c("coexist" = "#993399", # "#4477AA",
                 "yeast only" = "#FFCC33", # "#CCBB44",
                 "bacteria only" = "#333399", # "#AA3377",
                 "extinct" = "gray60")
outcome_shapes <- c("coexist" = 19,
                    "yeast only" = 19,
                    "bacteria only" = 19,
                    "extinct" = 4)

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "magenta")




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
                 g_yp = 0.0,  # 0.005, # pollinator-dependent for yeast (**)
                 g_b0 = 0.0,  # 0.02, # pollinator-independent for bacteria (**)
                 g_bp = 0.0,  # 0.002,  # pollinator-dependent for bacteria
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
                 dt = 0.1,
                 n_reps = 100,
                 season_len = 150,
                 max_t = 3000)


    # Arguments that aren't vectors:
    non_vecs <- c("max_t", "dt", "X", "u", "n_sigma", "n_reps",
                  "season_len", "season_surv", "season_sigma")

    for (n in names(args)[!names(args) %in% non_vecs]) {
        args[[n]] <- rep(args[[n]], np)
    }

    stopifnot(inherits(other_args, "list"))
    stopifnot(! is.null(names(other_args)) && ! any(names(other_args) == ""))

    for (n in names(other_args)) {
        if (! n %in% non_vecs && length(other_args[[n]]) == 1) {
            other_args[[n]] <- rep(other_args[[n]], np)
        }
        args[[n]] <- other_args[[n]]
    }


    args[["Y0"]] <- seq(0.1, 0.4, length.out = np)
    args[["B0"]] <- 0.5 - args[["Y0"]]
    stopifnot(args[["n_reps"]] %% 1 == 0)
    args[["n_reps"]] <- as.integer(args[["n_reps"]])

    if (no_error) {
        args[["n_sigma"]] <- NULL
    }

    run_df <- do.call(run_fun, args) |>
        as_tibble()

    if (is.null(args[["season_len"]])) {
        sl <- args[["max_t"]] * 0.1
    } else sl <- args[["season_len"]]
    if (is.null(args[["n_reps"]])) {
        nr <- 1L
    } else nr <- args[["n_reps"]]

    attr(run_df, "season_len") <- sl
    attr(run_df, "n_reps") <- nr
    attr(run_df, "dt") <- args[["dt"]]

    return(run_df)

}


summ_run <- function(even_run_output, .threshold = 1e-6) {

    sl <- attr(even_run_output, "season_len")
    nr <- attr(even_run_output, "n_reps")

    even_run_output |>
        # Take median through time (for each rep and plant) for last season:
        filter(t > max(t) - sl) |>
        group_by(rep, p) |>
        summarize(Y = median(Y),
                  B = median(B),
                  P = median(P),
                  .groups = "drop") |>
        # Take summed Y and B (across entire landscape) for each rep:
        group_by(rep) |>
        summarize(Y = sum(Y), B = sum(B)) |>
        mutate(outcome = case_when(Y > .threshold & B > .threshold ~ "coexist",
                                   Y > .threshold & B < .threshold ~ "yeast only",
                                   Y < .threshold & B > .threshold ~ "bacteria only",
                                   TRUE ~ "extinct") |>
                   factor(levels = names(outcome_pal))) |>
        split(~ outcome, drop = FALSE) |>
        imap_dfr(\(x, n) tibble(outcome = n, prob = nrow(x) / nr)) |>
        mutate(outcome = factor(outcome, levels = names(outcome_pal)))

}



# These are the bounds of where coexistence starts to occur, plus the middle
# of these bounds (when d_b0 = 0.3 and u = 1):
d_yp__ = c(`bacteria only` = 0.886308, coexist = 1.17755, `yeast only` = 1.46879) + c(-1, 0, 1) * 0.1




if (! file.exists(stoch_sims_file)) {

    # Takes ~30 min
    stoch_sim_df <- crossing(.u = 0:10,
                       .d_yp = d_yp__,
                       .n_sigma = c(100, 200, 400),
                       .season_surv = c(0.05, 0.1, 0.2),
                       .season_sigma = c(0, 10),
                       .np = c(2, 10)) |>
        # Don't do this in parallel bc landscape_constantF_stoch_ode is already
        # doing that
        pmap_dfr(\(.u, .d_yp, .n_sigma, .season_surv, .season_sigma, .np) {
            list(u = .u, d_yp = .d_yp, n_sigma = .n_sigma,
                 season_surv = .season_surv, season_sigma = .season_sigma) |>
                even_run(np = .np) |>
                summ_run() |>
                mutate(u = .u, d_yp = .d_yp, n_sigma = .n_sigma,
                       season_surv = .season_surv,
                       season_sigma = .season_sigma,
                       n_plants = .np)
        })
    write_rds(stoch_sim_df, stoch_sims_file)

} else {

    stoch_sim_df <- read_rds(stoch_sims_file)

}











season_sigma__ <- 0
n_plants__ <- 10

map(names(d_yp__), \(n) {
    p <- stoch_sim_df |>
        filter(d_yp == d_yp__[[n]],
               season_sigma == season_sigma__,
               n_plants == n_plants__) |>
        mutate(n_sigma = paste("n[sigma] ==", n_sigma),
               season_surv = paste("s ==", season_surv)) |>
        ggplot(aes(u, prob, color = outcome)) +
        ggtitle(paste("Outcome without stochasticity:\n", n)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
        geom_point(aes(shape = outcome)) +
        geom_line() +
        facet_grid(season_surv ~ n_sigma, labeller = label_parsed) +
        scale_y_continuous("Percent of simulations",
                           labels = scales::label_percent()) +
        scale_x_continuous("u", breaks = seq(0, 10, 2.5),
                           labels = c("0", "", "5", "", "10")) +
        scale_shape_manual(NULL, values = outcome_shapes) +
        scale_color_manual(NULL, values = outcome_pal) +
        theme(plot.title = element_text(size = 12,
                                        margin = margin(0,0,0,b=6)))
    if (n != names(d_yp__)[[1]]) {
        p <- p + theme(axis.title.y = element_blank(),
                       axis.text.y = element_blank())
    }
    if (n != tail(names(d_yp__), 1)) {
        p <- p + theme(strip.text.y = element_blank())
    }
    return(p)
}) |>
    c(list(guides = "collect")) |>
    do.call(what = wrap_plots)





# =============================================================================*
# =============================================================================*

# Open system sims ----

# =============================================================================*
# =============================================================================*



# Summarize for an open system where extinctions are temporary.
# In this case, output Bray-Curtis dissimilarity and Shannon diversity index
summ_open_run <- function(even_run_output) {

    .np <- length(unique(even_run_output$p))

    even_run_output |>
        group_by(rep) |>
        summarize(dive = diversity_vector(Y, B, group_size = .np),
                  diss = dissimilarity_vector(Y, B, group_size = .np, TRUE))

}





if (! file.exists(open_stoch_sims_file)) {

    # Takes ~70 min
    open_stoch_sim_df <- crossing(.u = 0:10,
                                  .d_yp = d_yp__,
                                  .n_sigma = c(100, 200, 400),
                                  .season_surv = c(0.05, 0.1, 0.2),
                                  .season_sigma = c(0, 10),
                                  .np = c(2, 10)) |>
        # Don't do this in parallel bc landscape_constantF_stoch_ode is already
        # doing that
        pmap_dfr(\(.u, .d_yp, .n_sigma, .season_surv, .season_sigma, .np) {
            list(u = .u, d_yp = .d_yp, n_sigma = .n_sigma,
                 season_surv = .season_surv, season_sigma = .season_sigma,
                 g_yp = 0.005, g_b0 = 0.02, g_bp = 0.002) |>
                even_run(np = .np) |>
                summ_open_run() |>
                mutate(u = .u, d_yp = .d_yp, n_sigma = .n_sigma,
                       season_surv = .season_surv,
                       season_sigma = .season_sigma,
                       n_plants = .np)
        })
    write_rds(open_stoch_sim_df, open_stoch_sims_file)

} else {

    open_stoch_sim_df <- read_rds(open_stoch_sims_file)

}



season_sigma__ <- 0
n_plants__ <- 2L




map(names(d_yp__), \(n) {
    dd <- open_stoch_sim_df |>
        filter(d_yp == d_yp__[[n]],
               season_sigma == season_sigma__,
               n_plants == n_plants__) |>
        mutate(n_sigma = paste("n[sigma] ==", n_sigma),
               season_surv = paste("s ==", season_surv)) |>
        pivot_longer(dive:diss) |>
        mutate(name = factor(name, levels = c("diss", "dive"),
                             labels = c("dissimilarity", "diversity")))
    dds <- dd |>
        group_by(u, d_yp, n_sigma, season_surv, season_sigma, name) |>
        summarize(value = median(value), .groups = "drop")
    p <- dd |>
        ggplot(aes(u, value, color = name)) +
        ggtitle(paste("Outcome without stochasticity:\n", n)) +
        geom_point(shape = 1, alpha = 0.1) +
        geom_line(data = dds, linewidth = 1) +
        facet_grid(season_surv ~ n_sigma, labeller = label_parsed) +
        scale_y_continuous("Measure value") +
        scale_x_continuous("u", breaks = seq(0, 10, 2.5),
                           labels = c("0", "", "5", "", "10")) +
        scale_color_viridis_d(NULL, option = "magma", begin = 0.3, end = 0.8) +
        theme(plot.title = element_text(size = 12,
                                        margin = margin(0,0,0,b=6)))
    if (n != names(d_yp__)[[1]]) {
        p <- p + theme(axis.title.y = element_blank(),
                       axis.text.y = element_blank())
    }
    if (n != tail(names(d_yp__), 1)) {
        p <- p + theme(strip.text.y = element_blank())
    }
    return(p)
}) |>
    c(list(guides = "collect")) |>
    do.call(what = wrap_plots)






# u__ <- 0
#
# z <- list(u = u__, d_yp = d_yp__[[2]], n_sigma = 200,
#           n_reps = 1,
#      g_yp = 0.005,
#      g_b0 = 0.02,
#      g_bp = 0.002,
#      season_surv = 0.1, season_sigma = 0) |>
#     even_run(np = 2L)
# zz <- z |>
#     group_by(t) |>
#     summarize(Y = sum(Y), B = sum(B))
#
#
# z |>
#     select(-P) |>
#     pivot_longer(Y:B, names_to = "species", values_to = "density") |>
#     mutate(species = factor(species, levels = c("Y", "B"),
#                             labels = c("yeast", "bacteria")),
#            p = factor(as.integer(p))) |>
#     ggplot(aes(t, density, color = species)) +
#     geom_line() +
#     facet_wrap(~ p) +
#     scale_color_manual(values = spp_pal)
#
# zz |>
#     pivot_longer(Y:B, names_to = "species", values_to = "density") |>
#     mutate(species = factor(species, levels = c("Y", "B"),
#                             labels = c("yeast", "bacteria"))) |>
#     ggplot(aes(t, density, color = species)) +
#     geom_line() +
#     scale_color_manual(values = spp_pal)




