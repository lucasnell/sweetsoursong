
library(sweetsoursong)
library(tidyverse)


b_curve <- function(s0, h, ...) {
    curve(s0^h / (s0^h + x^h), 0, 1,
          xlab =  "Proportion bacteria",
          ylab = "Pollinator density at plant",
          ylim = c(0, 1), lwd = 2,
          ...)
}
f_curve <- function(f0, u, ...) {
    curve(x^u / (f0^u + x^u), 0, 100,
          xlab =  "Total flowers",
          ylab = "Pollinator density at plant",
          ylim = c(0, 1), lwd = 2,
          ...)
}


high_low_run <- function(other_args) {

    if (inherits(other_args, "data.frame")) {
        stopifnot(nrow(other_args) == 1L)
        other_args <- as.list(other_args)
    }

    args <- list(dt = 0.1,
                 max_t = 100,
                 m = 0.1,
                 R = 10)

    stopifnot(inherits(other_args, "list"))
    stopifnot(! is.null(names(other_args)) && ! any(names(other_args) == ""))

    if ("max_t" %in% names(other_args)) {
        args[["max_t"]] <- other_args[["max_t"]]
        other_args[["max_t"]] <- NULL
    }

    stopifnot(! any(names(other_args) %in% names(args)))

    all_runs_args <-list("low Y" = list(Y0 = 1,  B0 = 10, N0 = 89),
                         "low B" = list(Y0 = 10, B0 = 1,  N0 = 89),
                         "both low" = list(Y0 = 1, B0 = 1, N0 = 98),
                         "both high" = list(Y0 = 10, B0 = 10, N0 = 80))

    all_runs <- all_runs_args |>
        imap_dfr(\(x, n) {
            c(args, other_args, x) |>
                do.call(what = one_plant_ode) |>
                as_tibble() |>
                mutate(run = n)
        }) |>
        select(run, everything()) |>
        mutate(run = factor(run, levels = names(all_runs_args)))

    return(all_runs)

}



tibble(# ----------
    # within-plant dispersal rates:
    d_yp = 1.5,  # pollinator-dependent for yeast (**)
    d_b0 = 0.3,  # pollinator-independent for bacteria (**)
    d_bp = 0.4,  # pollinator-dependent for bacteria (**)
    # -------------
    # rates of immigration from non-focal-plant sources:
    g_yp = 0.005, # pollinator-dependent for yeast (**)
    g_b0 = 0.02, # pollinator-independent for bacteria (**)
    g_bp = 0.001,  # pollinator-dependent for bacteria
    # -------------
    # half saturation ratios:
    L_0 = 0.01,  # for P -> dispersal (**)
    s_0 = 0.5,   # for B/F -> P (??)
    f_0 = 0.5,   # for F -> P (??)
    # -------------
    # shape exponents:
    h = 3,  # for B/F -> P (??)
    u = 2,  # for F -> P (??)
    # -------------
    # others:
    P_max = 1.5,  # maximum possible pollinator density
    q = 0.5,      # relative strength of B/F -> P versus F -> P (??)
    F_tilde = 1e3,# number of nearby flowers other than focal plant (??)
    #
    # (**)  = value from Song et al. (submitted)
    # (??) = value should be varied bc I have no clue what to use
    #
    max_t = 250
    ) |>
    high_low_run() |>
    (\(x) {
        print(group_by(x, run) |>
                  summarize(Y = tail(Y, 1),
                      B = tail(B, 1),
                      N = tail(N, 1),
                      P = tail(P, 1)))
        return(x)
    })() |>
    mutate(P = P * 50) |>
    pivot_longer(Y:P, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B", "N", "P"),
                         labels = c("yeast", "bacteria", "non-colonized",
                                    "pollinators"))) |>
    ggplot(aes(t, density)) +
    geom_hline(yintercept = 0, linewidth = 1, linetype = "22", color = "gray70") +
    geom_line(aes(color = type), linewidth = 1) +
    facet_wrap(~ run, ncol = 1) +
    xlab("Time (days)") +
    scale_y_continuous("flower-type density",
                       sec.axis = sec_axis(~ . / 50, "pollinator density")) +
    scale_color_viridis_d(NULL, begin = 0.1, option = "H")
