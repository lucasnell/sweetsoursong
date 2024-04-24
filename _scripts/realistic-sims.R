
library(MASS)  # mvrnorm
library(terra) # autocor
library(spatstat.random) # rMatClust
library(sweetsoursong)
library(tidyverse)
library(patchwork)
library(gganimate)


spp_pal <- c("non-colonized" = "gray60",
             "yeast" = "#FFCC33",
             "bacteria" = "#333399",
             "pollinators" = "firebrick3")



# Flowering through time fit to normal distributions for 30 plants in JRBP:
norm_phen_fits <- read_rds("_data/norm-phen-fits.rds") |>
    split(~ plant) |>
    map_dfr(\(x) {
        stopifnot(nrow(x) == 1)
        t <- 0:365
        r <- x$sum * dnorm(t, x$mean, x$sd)
        return(mutate(x, min_t = min(t[r >= 1]), max_t = max(t[r >= 1])))
    })

# Bounds on a flowering season:
(flower_start <- norm_phen_fits$min_t |> min())
(flower_stop <- norm_phen_fits$max_t |> max())



get_phenology <- function(.n_plants) {
    if (.n_plants > nrow(norm_phen_fits)) {
        vcvm <- matrix(0, 3, 3)
        vars <- list("mean" = log(norm_phen_fits[["mean"]]),
                     "sd" = log(norm_phen_fits[["sd"]]),
                     "sum" = log(norm_phen_fits[["sum"]]))
        for (i in 1:3) vcvm[i,i] <- var(vars[[i]])
        for (i in 2:3) {
            for (j in 1:(i-1)) {
                covij <- cov(vars[[i]], vars[[j]])
                vcvm[i,j] <- vcvm[j,i] <- covij
            }
        }
        means <- sapply(vars, mean)
        z <- mvrnorm(.n_plants, mu = means, Sigma = vcvm)
        # make sure they don't exceed observed bounds:
        for (n in names(vars)) {
            z[z[,n] > max(vars[[n]]), n] <- max(vars[[n]])
            z[z[,n] < min(vars[[n]]), n] <- min(vars[[n]])
            z[,n] <- exp(z[,n])
        }
        args <- list(R_hat = z[,"sum"], sigma = z[,"sd"],
                     mu = z[,"mean"] - flower_start)
    } else {
        plt_idx <- sample.int(nrow(norm_phen_fits), .n_plants,
                              replace = FALSE)
        args <- list(R_hat = norm_phen_fits[["sum"]][plt_idx],
                     sigma = norm_phen_fits[["sd"]][plt_idx],
                     mu = norm_phen_fits[["mean"]][plt_idx] - flower_start)
    }
    return(args)
}





one_run <- function(...) {

    other_args <- list(...)

    if (inherits(other_args, "data.frame")) {
        stopifnot(nrow(other_args) == 1L)
        other_args <- as.list(other_args)
    }
    stopifnot(inherits(other_args, "list"))
    stopifnot(!is.null(names(other_args)) && all(names(other_args) != ""))
    stopifnot(all(names(other_args) %in% names(formals(landscape_season_ode))))
    stopifnot("z" %in% names(other_args))
    if (is.list(other_args$z)) other_args$z <- other_args$z[[1]]
    stopifnot("Y0" %in% names(other_args))
    stopifnot("B0" %in% names(other_args))

    stopifnot(!is.null(nrow(other_args$z)) && nrow(other_args$z) >= 2)
    .n_plants <- nrow(other_args$z)

    for (n in names(other_args)) {
        x <- other_args[[n]]
        if (n %in% c("S_0", "q", "w", "add_F")) {
            if (length(x) != 1) {
                stop(sprintf("\nERROR: Argument %s must be length %i, but it's %i.\n",
                             n, 1, length(x)))
            }
        } else if (n == "z") {
            if (! identical(dim(x), rep(.n_plants, 2)) || ! all(x >= 0) ) {
                stop("Argument z must be a square numeric matrix with all values >= 0")
            }
        } else {
            if (length(x) == 1) {
                other_args[[n]] <- rep(x, .n_plants)
            } else if (length(x) != .n_plants) {
                stop(sprintf("\nERROR: Argument %s must be length %i, but it's %i.\n",
                             n, .n_plants, length(x)))
            }
        }
        if (! is.numeric(x)) stop(paste("\nERROR: Argument", n, "must be numeric.\n"))
    }

    arg_list <- list(# ----------*
                     # within-plant dispersal rates:
                     d_yp = rep(1.5, .n_plants),  # pollinator-dependent for Y (**)
                     d_b0 = rep(0.3, .n_plants),  # pollinator-independent for B (**)
                     d_bp = rep(0.4, .n_plants),  # pollinator-dependent for B (**)
                     # -------------*
                     # rates of immigration from non-focal-plant sources:
                     g_yp = rep(0.005, .n_plants), # pollinator-dependent for Y (**)
                     g_b0 = rep(0.02, .n_plants), # pollinator-independent for B (**)
                     g_bp = rep(0.001, .n_plants),  # pollinator-dependent for B
                     # -------------*
                     # half saturation ratios:
                     L_0 = rep(0.01, .n_plants),  # for P/F -> dispersal (**)
                     S_0 = 0.6,   # for B/F -> P (??)
                     # -------------*
                     # others:
                     P_max = rep(.n_plants / 2, .n_plants), # maximum P density
                     X = rep(0, .n_plants), # attraction to non-focal flowers (??)
                     q = 0, # affects the strength of F -> P (??)
                     w = 0, # how quickly effects of nearby focal flowers affect P (??)
                     add_F = 1, # F at which to add Y0 and B0
                     #
                     # (**)  = value from Song et al. (submitted)
                     # (??) = value should be varied bc I have no clue what to use
                     #
                     dt = 0.1,
                     max_t = flower_stop - flower_start + 20,
                     m = rep(0.2, .n_plants))

    arg_list <- c(arg_list, get_phenology(.n_plants))
    for (n in names(other_args)) arg_list[[n]] <- other_args[[n]]

    sim_df <- do.call(landscape_season_ode, arg_list) |>
                as_tibble() |>
        mutate(p = factor(as.integer(p) + 1L))

    return(sim_df)

}



time_series <- function(run_df, main = NULL, ylim_F = NULL, ylim_P = NULL) {
    dd <- run_df |>
        pivot_longer(Y:P, names_to = "type", values_to = "density") |>
        mutate(type = factor(type, levels = c("Y", "B", "N", "P"),
                             labels = c("yeast", "bacteria", "non-colonized",
                                        "pollinators")),
               id = interaction(p, type, drop = TRUE)) |>
        filter(t %% 1 == 0)
    if (length(unique(run_df$p)) > 2) {
        .a <- 0.25
        if (length(unique(run_df$p)) > 50) .a <- 0.1
        dd <- dd |>
            filter(type != "non-colonized")
    } else {
        .a <- 1
    }

    p1 <- dd |>
        filter(type != "pollinators") |>
        ggplot(aes(t, density)) +
        geom_hline(yintercept = 0, linewidth = 1, linetype = "22", color = "gray70") +
        geom_line(aes(color = type, group = id), linewidth = 1, alpha = .a) +
        xlab("Time (days)") +
        ggtitle(main) +
        scale_y_continuous("flower-type density", limits = ylim_F) +
        scale_color_manual(NULL, values = spp_pal, guide = "none") +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5))
    p2 <- dd |>
        filter(type == "pollinators") |>
        ggplot(aes(t, log1p(density))) +
        geom_hline(yintercept = 0, linewidth = 1, linetype = "22", color = "gray70") +
        geom_line(aes(group = id), linewidth = 1, alpha = .a,
                  color = spp_pal[["pollinators"]]) +
        xlab("Time (days)") +
        scale_y_continuous("log(1 + pollinators)", limits = ylim_P)
    p1 / p2
}


Y_hist <- function(runs_df, xlim = NULL) {
    runs_df |>
        filter((Y + B + N) > 0) |>
        group_by(q, p) |>
        summarize(Y = mean(100 * Y / (Y + B + N)), .groups = "drop") |>
        ggplot(aes(Y)) +
        geom_histogram(bins = 20, fill = spp_pal[["yeast"]]) +
        # scale_x_continuous("mean % yeast") +
        # coord_cartesian(xlim = xlim_Y)
        scale_x_continuous("mean % yeast", limits = xlim) +
        facet_wrap(~ q, ncol = 1, scales = "free_y",
                   labeller = \(x) label_both(x, sep = " = "))
}

P_hist <- function(run_df, xlim = NULL) {
    run_df |>
        group_by(q, p) |>
        summarize(P = mean(log1p(P)), .groups = "drop") |>
        ggplot(aes(P)) +
        geom_histogram(bins = 20, fill = spp_pal[["pollinators"]]) +
        scale_x_continuous("mean log(1 + pollinators)", limits = xlim) +
        facet_wrap(~ q, ncol = 1, scales = "free_y",
                   labeller = \(x) label_both(x, sep = " = "))
}



animated_landscape <- function(run_df, xy_df, filename = NULL) {

    bounds <- list(x = range(xy_df$x) + c(-1,1) * 0.1 * diff(range(xy_df$x)),
                   y = range(xy_df$y) + c(-1,1) * 0.1 * diff(range(xy_df$y)))

    p <- run_df |>
        filter(t %% 1 == 0) |>
        mutate(t = as.integer(t),
               x = map_dbl(p, \(pp) xy_df$x[[as.integer(paste(pp))]]),
               y = map_dbl(p, \(pp) xy_df$y[[as.integer(paste(pp))]]),
               nf = Y + B + N,
               Y = 100 * Y / nf) |>
        ggplot(aes(x, y)) +
        geom_point(aes(size = nf, color = Y)) +
        scale_color_viridis_c("% yeast", limits = c(0, 100)) +
        scale_size("# flowers", limits = c(0, NA)) +
        coord_equal(xlim = bounds$x, ylim = bounds$y)

    ap <- p + transition_time(t) +
        labs(title = "Day: {frame_time}")

    if (!is.null(filename)) {
        anim_save(sprintf("_figures/%s.gif", filename), ap,
                  height = 4, width = 6, units = "in", res = 150)
        invisible(NULL)
    }

    return(ap)
}



#' ========================================================================
#' ========================================================================
# RANDOM LANDSCAPE ----
#' ========================================================================
#' ========================================================================

np <- 100

set.seed(678902345)
xy_df <- tibble(x = runif(np), y = runif(np))
dist_mat <- make_dist_mat(xy_df)
Y0 <- runif(100, 0, 0.5)
B0 <- 0.5 - Y0
phen <- get_phenology(np)

# xy_df |>
#     ggplot(aes(x, y)) +
#     geom_point() +
#     theme(axis.title = element_blank(),
#           axis.text = element_blank(),
#           axis.ticks = element_blank()) +
#     coord_equal(xlim = c(0, 1), ylim = c(0, 1))

pmap_dfr(c(phen, list(p = 1:np)), \(R_hat, sigma, mu, p) {
    tibble(plant = factor(p, levels = 1:np),
           time = 0:(flower_stop - flower_start + 20),
           r = R_hat * dnorm(time, mu, sigma))
}) |>
    ggplot(aes(time, r, group = plant)) +
    geom_line(alpha = 0.25, linewidth = 1) +
    theme(axis.title = element_blank())

# showing effects of q:
crossing(x = 1:100, h = c(0.5, 1:3)) |>
    mutate(sum_h = map_dbl(h, \(.h) sum((1:100)^.h)), z = x^h / sum_h) |>
    ggplot(aes(x, z, color = factor(h))) +
    geom_line(linewidth = 1) +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_brewer(palette = "Dark2", guide = "none")


runs_df <- map_dfr(c(0.5, 1, 2, 3),
                   \(.q) {
                       one_run(Y0 = Y0 * 0, B0 = B0 * 0,
                               z = dist_mat,
                               X = 25^(1/.q),
                               q = .q,
                               w = 100,
                               R_hat = phen$R_hat,
                               sigma = phen$sigma,
                               add_F = 1,
                               mu = phen$mu) |>
                           mutate(q = factor(.q, levels = c(0.5, 1, 2, 3))) |>
                           dplyr::select(q, everything())
                   })
ts_ps <- runs_df |>
    split(~ q) |>
    map(\(x) {
        p <- time_series(x, ylim_F = c(0, 250), ylim_P = c(0, 4),
                         main = sprintf("q = %s", paste(x$q[[1]])))
        return(p)
    })

do.call(wrap_plots, c(ts_ps, list(ncol = 2)))



Y_p <- runs_df |> Y_hist()
P_p <- runs_df |> P_hist()



#' ========================================================================
#' ========================================================================
# CLUSTERED LANDSCAPE ----
#' ========================================================================
#' ========================================================================



cluster_radii <- c(0.02, 0.1, 1)
w_vals <- c(1, 10, 100)

set.seed(13456978)
clust_sims <- map(cluster_radii,  \(.s) rMatClust(kappa = 25, scale = .s, mu = 10, nsim = 100))
for (i in 1:length(clust_sims)) stopifnot(all(map_int(clust_sims[[i]], \(x) x$n) >= np))
# force all to have exactly `np` plants:
for (i in 1:length(clust_sims)) {
    for (j in 1:length(clust_sims[[i]])) {
        idx <- sample.int(clust_sims[[i]][[j]][["n"]], np, replace = FALSE)
        clust_sims[[i]][[j]][["x"]] <- clust_sims[[i]][[j]][["x"]][idx]
        clust_sims[[i]][[j]][["y"]] <- clust_sims[[i]][[j]][["y"]][idx]
        clust_sims[[i]][[j]][["n"]] <- np
    }
}
for (i in 1:length(clust_sims)) stopifnot(all(map_int(clust_sims[[i]], \(x) x$n) == np))
names(clust_sims) <- cluster_radii


cl_runs_df <- map2_dfr(rep(names(clust_sims), length(w_vals)),
                       rep(w_vals, each = length(cluster_radii)),
                       \(.s, .w) {
                           cs <- clust_sims[[.s]][[1]]
                           dist_mat__ <- make_dist_mat(cs$x, cs$y)
                           one_run(Y0 = Y0 * 0, B0 = B0 * 0,
                                   z = dist_mat__,
                                   X = 25,
                                   q = 1,
                                   w = .w,
                                   R_hat = phen$R_hat,
                                   sigma = phen$sigma,
                                   mu = phen$mu) |>
                               mutate(s = factor(.s, levels = names(clust_sims)),
                                      w = factor(.w, levels = w_vals)) |>
                               dplyr::select(s, w, everything())
                       })



cl_ts_ps <- cl_runs_df |>
    filter(s == 1) |>
    split(~ w) |>
    map(\(x) time_series(x, ylim_F = c(0, 200), ylim_P = c(0, 12),
                         main = sprintf("w = %s", paste(x$w[[1]]))))


do.call(wrap_plots, c(cl_ts_ps, list(ncol = 1)))

cl_runs_df |> Y_hist2()
cl_runs_df |> P_hist2()


Y_hist2 <- function(cl_runs_df, xlim = NULL) {
    cl_runs_df |>
        filter((Y + B + N) > 0, s == 1) |>
        group_by(s, w, p) |>
        summarize(Y = mean(100 * Y / (Y + B + N)), .groups = "drop") |>
        ggplot(aes(Y)) +
        geom_histogram(bins = 20, fill = spp_pal[["yeast"]]) +
        scale_x_continuous("mean % yeast", limits = xlim) +
        facet_wrap(~ w, ncol = 1, labeller = \(x) label_both(x, sep = " = "))
}

P_hist2 <- function(cl_runs_df, xlim = NULL) {
    cl_runs_df |>
        filter(s == 1) |>
        group_by(s, w, p) |>
        summarize(P = mean(P), .groups = "drop") |>
        ggplot(aes(P)) +
        geom_histogram(bins = 20, fill = spp_pal[["pollinators"]]) +
        scale_x_continuous("mean pollinator density", limits = xlim) +
        facet_wrap(~ w, ncol = 1, labeller = \(x) label_both(x, sep = " = "))
}












# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# OLD CODE BELOW
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



high_low_run <- function(other_args, .n_plants = 30L) {

    if (inherits(other_args, "data.frame")) {
        stopifnot(nrow(other_args) == 1L)
        other_args <- as.list(other_args)
    }
    stopifnot(!is.null(names(other_args)) && all(names(other_args) != ""))
    stopifnot(all(names(other_args) %in% names(formals(landscape_season_ode))))

    stopifnot(is.numeric(.n_plants) && length(.n_plants) == 1 && .n_plants %% 1 == 0)
    .n_plants <- as.integer(.n_plants)
    if (.n_plants < 2L) {
        stop("\nERROR: Must be at least 2 plants.\n")
    }

    if (.n_plants <= nrow(norm_phen_fits)) {
        plt_idx <- 1:.n_plants
    } else plt_idx <- sample.int(nrow(norm_phen_fits), .n_plants, replace = TRUE)

    args <- list(dt = 0.1,
                 max_t = max(norm_phen_fits$max_t) - min(norm_phen_fits$min_t) + 20,
                 m = rep(0.2, .n_plants),
                 R_hat = norm_phen_fits[["sum"]][plt_idx],
                 sigma = norm_phen_fits[["sd"]][plt_idx],
                 mu = norm_phen_fits[["mean"]][plt_idx] - min(norm_phen_fits$min_t))

    stopifnot(inherits(other_args, "list"))
    stopifnot(! is.null(names(other_args)) && ! any(names(other_args) == ""))

    if ("max_t" %in% names(other_args)) {
        args[["max_t"]] <- other_args[["max_t"]]
        other_args[["max_t"]] <- NULL
    }

    for (n in names(other_args)) {
        x <- other_args[[n]]
        if (n %in% c("S_0", "q", "w", "add_F")) {
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

    hi_seq <- seq(0.4,  0.25, length.out = .n_plants)
    lo_seq <- seq(0.25, 0.1,  length.out = .n_plants)
    eq_seq <- seq(0.4, 0.1,   length.out = .n_plants)

    all_runs_args <- list("high Y" = list(Y0 = hi_seq, B0 = rev(lo_seq)),
                          "high B" = list(Y0 = lo_seq, B0 = rev(hi_seq)),
                          "equal" =  list(Y0 = eq_seq, B0 = rev(eq_seq)))

    all_runs <- all_runs_args |>
        imap_dfr(\(x, n) {
            c(args, other_args, x) |>
                do.call(what = landscape_season_ode) |>
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
    })

# xy_df |>
#     ggplot(aes(x, y)) +
#     geom_point() +
#     coord_equal(xlim = c(0, 10), ylim = c(0, 10))


# np <- 2L
dist_mat <- matrix(1, np, np) - diag(np)

for (i in 2:np) {
    for (j in 1:(i-1)) {
        dist_mat[i,j] <- dist_mat[j,i] <- sqrt((xy_df$x[i] - xy_df$x[j])^2 +
                                                   (xy_df$y[i] - xy_df$y[j])^2)
    }
}


make_plot <- function(hl_df) {
    max_YBN <- max(pmax(hl_df$Y, hl_df$B, hl_df$N))
    max_P <- max(hl_df$P)
    dd <- hl_df |>
        mutate(P = P * max_YBN / max_P) |>
        pivot_longer(Y:P, names_to = "type", values_to = "density") |>
        mutate(type = factor(type, levels = c("Y", "B", "N", "P"),
                             labels = c("yeast", "bacteria", "non-colonized",
                                        "pollinators")),
               id = interaction(run, p, type, drop = TRUE))
    if (length(unique(hl_df$p)) > 2) {
        .facet <- facet_wrap(~ run, ncol = 1)
        .a <- 0.25
        dd <- dd |>
            filter(t %% 1 == 0, type != "non-colonized")
    } else {
        .facet <- facet_grid(run ~ p)
        .a <- 1
    }

    dd |>
        ggplot(aes(t, density)) +
        geom_hline(yintercept = 0, linewidth = 1, linetype = "22", color = "gray70") +
        geom_line(aes(color = type, group = id), linewidth = 1, alpha = .a) +
        .facet +
        xlab("Time (days)") +
        scale_y_continuous("flower-type density",
                           sec.axis = sec_axis(~ . * max_P / max_YBN,
                                               "pollinator density")) +
        scale_color_manual(NULL, values = spp_pal, guide = "none")
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
       add_F = 1,          # when to add Y0 and B0
       #
       # (**)  = value from Song et al. (submitted)
       # (??) = value should be varied bc I have no clue what to use
       #
       )


four_args <- function(.S_0, .X, .q, .w) {
    hl_df <- arg_df |>
        mutate(S_0 = .S_0,
               X = .X,
               q = .q,
               w = .w) |>
        high_low_run(.n_plants = np)
    make_plot(hl_df) + ggtitle(parse(text = sprintf("S[0] == %s * ',' ~ X == %s * ',' ~ q == %s * ',' ~ w == %s", .S_0, .X, .q, .w)))
    # return(hl_df)
}

four_args(.S_0 = 0.6, .X = 25, .q = 0, .w = 10e3)
four_args(.S_0 = 0.6, .X = 25, .q = 1, .w = 10e3)
four_args(.S_0 = 0.6, .X = 25, .q = 2, .w = 10e3)

four_args(.S_0 = 0.6, .X = 25, .q = 2, .w = 10e3)
four_args(.S_0 = 0.6, .X = 25, .q = 2, .w = 1)
four_args(.S_0 = 0.6, .X = 25, .q = 2, .w = 0.1)





# x <- rbind(c( 6.33e-12, 2.59e-11, 0.0000000243),
#            c(7.37e-13, 3.01e-12, 0.00000000284))
# sweetsoursong:::landscape_weights(x = x,
#                   S_0 = 0.6,
#                   q = 2,
#                   X = rep(0.1 / np, 2),
#                   w = 10e3,
#                   z = matrix(1, 2, 2) - diag(2))



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
