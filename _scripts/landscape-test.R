
library(sweetsoursong)
library(tidyverse)
library(patchwork)


spp_pal <- c("non-colonized" = "gray60",
                 "yeast" = "#FFCC33",
                 "bacteria" = "#333399",
                 "pollinators" = "firebrick3")

high_low_run <- function(other_args, .n_plants = 2L) {

    if (inherits(other_args, "data.frame")) {
        stopifnot(nrow(other_args) == 1L)
        other_args <- as.list(other_args)
    }

    stopifnot(is.numeric(.n_plants) && length(.n_plants) == 1 && .n_plants %% 1 == 0)
    .n_plants <- as.integer(.n_plants)
    if (.n_plants < 2L) {
        stop("\nERROR: Must be at least 2 plants.\n")
    }

    args <- list(dt = 0.1,
                 max_t = 100,
                 m = rep(0.1, .n_plants),
                 R = rep(10, .n_plants),
                 N0 = rep(89, .n_plants))

    stopifnot(inherits(other_args, "list"))
    stopifnot(! is.null(names(other_args)) && ! any(names(other_args) == ""))

    if ("max_t" %in% names(other_args)) {
        args[["max_t"]] <- other_args[["max_t"]]
        other_args[["max_t"]] <- NULL
    }

    for (n in names(other_args)) {
        x <- other_args[[n]]
        if (n %in% c("S_0", "w", "q")) {
            if (length(x) != 1) {
                stop(sprintf("\nERROR: Argument %s must be length %i, but it's %i.\n",
                             n, 1, length(x)))
            }
        } else if (n == "z") {
            if (is.list(x)) {
                x <- x[[1]]
                other_args[[n]] <- x
            }
            if (! identical(dim(x), rep(.n_plants, 2)) || ! all(x >= 0) ) {
                stop("Argument z must be a square numeric matrix with all values >= 0")
            }
        } else {
            if (length(x) == 1) {
                x <- rep(x, .n_plants)
                other_args[[n]] <- x
            } else if (length(x) != .n_plants) {
                stop(sprintf("\nERROR: Argument %s must be length %i, but it's %i.\n",
                             n, .n_plants, length(x)))
            }
        }
        if (! is.numeric(x)) stop(paste("\nERROR: Argument", n, "must be numeric.\n"))
    }


    stopifnot(! any(names(other_args) %in% names(args)))

    all_runs_args <- list("high Y" =   list(Y0 = seq(10, 6, length.out = .n_plants),
                                            B0 = seq(1, 5, length.out = .n_plants)),
                          "high B" =   list(Y0 = seq(5, 1, length.out = .n_plants),
                                            B0 = seq(6, 10, length.out = .n_plants)),
                          "equal" = list(Y0 = seq(10, 1, length.out = .n_plants),
                                         B0 = seq(1, 10, length.out = .n_plants)))

    all_runs <- all_runs_args |>
        imap_dfr(\(x, n) {
            c(args, other_args, x) |>
                do.call(what = landscape_ode) |>
                as_tibble() |>
                mutate(run = n)
        }) |>
        select(run, everything()) |>
        mutate(run = factor(run, levels = names(all_runs_args)),
               p = factor(as.integer(p) + 1L))

    return(all_runs)

}


n_plt_per_patch <- 1:5 * 2L
patch_radius <- 0.5

np <- sum(n_plt_per_patch)


patch_centers <- tibble(x = c(1, 1, 5, 9, 9),
                        y = c(1, 9, 5, 1, 9))

xy_df <- map_dfr(1:nrow(patch_centers),
                 \(i) {
                     .x <- patch_centers[["x"]][[i]]
                     .y <- patch_centers[["y"]][[i]]
                     .np <- n_plt_per_patch[[i]]
                     .angles <- seq(0, 2, length.out = .np+1L)[-1] * pi
        tibble(x = .x + patch_radius * cos(.angles),
               y = .y + patch_radius * sin(.angles))
    }) |>
    #' Shuffling rows to make sure positions and proportion Y or B at start
    #' aren't correlated:
    slice_sample(prop = 1)

# xy_df |>
#     ggplot(aes(x, y)) +
#     geom_point() +
#     coord_equal(xlim = c(0, 10), ylim = c(0, 10))



dist_mat <- matrix(0, np, np)

for (i in 2:np) {
    for (j in 1:(i-1)) {
        dist_mat[i,j] <- dist_mat[j,i] <- sqrt((xy_df$x[i] - xy_df$x[j])^2 +
                                                   (xy_df$y[i] - xy_df$y[j])^2)
    }
}


make_plots <- function(hl_df) {
    hl_ts_p <- hl_df |>
        # (\(x) {
        #     print(group_by(x, run) |>
        #               filter(t == max(t)) |>
        #               summarize(Y = mean(Y),
        #                         B = mean(B),
        #                         N = mean(N),
        #                         Pmin = min(P),
        #                         Pmax = max(P),
        #                         .groups = "drop"))
        #
        #     return(x)
        # })() |>
        filter(t %% 1 == 0) |>
        mutate(P = P * 50 / max(hl_df$P)) |>
        pivot_longer(Y:P, names_to = "type", values_to = "density") |>
        mutate(type = factor(type, levels = c("Y", "B", "N", "P"),
                             labels = c("yeast", "bacteria", "non-colonized",
                                        "pollinators")),
               id = interaction(run, p, type, drop = TRUE)) |>
        filter(type != "non-colonized") |>
        ggplot(aes(t, density)) +
        geom_hline(yintercept = 0, linewidth = 1, linetype = "22", color = "gray70") +
        geom_line(aes(color = type, group = id), linewidth = 1, alpha = 0.25) +
        facet_wrap(~ run, ncol = 1) +
        xlab("Time (days)") +
        scale_y_continuous("flower-type density") + # ,
        # sec.axis = sec_axis(~ . * max(hl_df$P) / 100, "pollinator density")) +
        scale_color_manual(NULL, values = spp_pal, guide = "none")


    hl_fd_p <- hl_df |>
        filter(t == max(t)) |>
        mutate(P = P * 100 / max(hl_df$P)) |>
        select(run, p, Y, B, N, P) |>
        pivot_longer(Y:P, names_to = "type", values_to = "density") |>
        mutate(type = factor(type, levels = c("Y", "B", "N", "P"),
                             labels = c("yeast", "bacteria", "non-colonized",
                                        "pollinators"))) |>
        filter(type != "non-colonized") |>
        ggplot(aes(type, density)) +
        geom_hline(yintercept = 0, linewidth = 1, linetype = "22", color = "gray70") +
        geom_jitter(aes(color = type), height = 0, width = 0.2, alpha = 1,
                    shape = 1, size = 2) +
        facet_wrap( ~ run, ncol = 1) +
        scale_color_manual(values = spp_pal, guide = "none") +
        theme(axis.title.x = element_blank())

    hl_ts_p + hl_fd_p

}

arg_df <- tibble(# ----------
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
       L_0 = 0.01,  # for P/F -> dispersal (**)
       S_0 = 0.6,     # for B/F -> P (??)
       # -------------
       # others:
       P_max = np / 2,    # maximum possible pollinator density
       X = 0, # 0.01 / np,        # attraction to nearby flowers other than focal plant (??)
       q = 0,            # affects the strength of F -> P (??)
       w = 0,        # how quickly effects of nearby focal flowers affect pollinators (??)
       z = list(dist_mat), # distance matrix
       #
       # (**)  = value from Song et al. (submitted)
       # (??) = value should be varied bc I have no clue what to use
       #
       max_t = 100)


four_args <- function(.S_0, .X, .q, .w) {
    hl_df <- arg_df |>
        mutate(S_0 = .S_0,
               X = .X,
               q = .q,
               w = .w) |>
        high_low_run(.n_plants = np)
    make_plots(hl_df) + plot_annotation(title = parse(text = sprintf("S[0] == %s * ',' ~ X == %s * ',' ~ q == %s * ',' ~ w == %s", .S_0, .X, .q, .w)))
}

four_args(.S_0 = 0.6, .X = 0, .q = 0, .w = 0)
four_args(.S_0 = 10e3, .X = 0, .q = 0, .w = 0)
#
four_args(.S_0 = 0.6, .X = 0, .q = 0, .w = 0)
four_args(.S_0 = 0.6, .X = 0.001, .q = 0, .w = 0)
#
four_args(.S_0 = 0.6, .X = 0, .q = 0, .w = 0)
four_args(.S_0 = 20, .X = 0, .q = 0, .w = 0)







hl_df |>
    filter(t == max(t)) |>
    select(run, p, Y, B, N, P) |>
    mutate(x = map_dbl(p, \(.p) xy_df$x[as.integer(paste(.p))]),
           y = map_dbl(p, \(.p) xy_df$y[as.integer(paste(.p))])) |>
    ggplot(aes(x, y)) +
    geom_point(aes(color = log(Y))) +
    facet_wrap(~ run) +
    coord_equal()


hl_df |>
    filter(t == max(t)) |>
    select(run, p, Y, B, N, P) |>
    mutate(Y0 = map2_dbl(p, run, \(.p, .r) {
        pp <- as.integer(paste(.p))
        all_runs_args <- list("high Y" =   list(Y0 = seq(10, 6, length.out = np),
                                                B0 = seq(1, 5, length.out = np)),
                              "high B" =   list(Y0 = seq(5, 1, length.out = np),
                                                B0 = seq(6, 10, length.out = np)),
                              "equal" = list(Y0 = seq(10, 1, length.out = np),
                                             B0 = seq(1, 10, length.out = np)))
        return(all_runs_args[[.r]][["Y0"]][[pp]])
    })) |>
    ggplot(aes(Y0, Y)) +
    geom_point()
