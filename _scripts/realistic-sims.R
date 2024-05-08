
suppressPackageStartupMessages({
    library(MASS)  # mvrnorm
    library(terra) # autocor
    library(spatstat.random) # rMatClust
    library(sweetsoursong)
    library(tidyverse)
    library(patchwork)
    library(gganimate)
    library(parallel)
})



spp_pal <- c("non-colonized" = "gray60",
             "yeast" = "#FFCC33",
             "bacteria" = "#333399",
             "pollinators" = "firebrick3",
             "total" = "black")


if (file.exists(".Rprofile")) source(".Rprofile")

options("mc.cores" = max(1L, detectCores() - 2L))


# Flowering through time fit to normal distributions for 30 plants in JRBP:
phen_fits <- read_rds("_data/norm-phen-fits.rds") |>
    split(~ plant) |>
    map_dfr(\(x) {
        stopifnot(nrow(x) == 1)
        t <- 0:365
        r <- x$sum * dnorm(t, x$mean, x$sd)
        return(mutate(x, min_t = min(t[r >= 1]), max_t = max(t[r >= 1])))
    })

# Bounds on a flowering season:
flower_start <- phen_fits$min_t |> min()
flower_stop <- phen_fits$max_t |> max()




get_phenology <- function(.n_plants) {
    stopifnot(is.numeric(.n_plants) && all(.n_plants %% 1 == 0))
    # .n_plants = 100L
    # rm(.n_plants, vcvm, vars, i, j, covij, means, z, n, args, plt_idx)
    if (.n_plants > nrow(phen_fits)) {
        vcvm <- matrix(0, 3, 3)
        vars <- list("mean" = log(phen_fits[["mean"]]),
                     "sd" = log(phen_fits[["sd"]]),
                     "sum" = log(phen_fits[["sum"]]))
        for (i in 1:3) vcvm[i,i] <- var(vars[[i]])
        for (i in 2:3) {
            for (j in 1:(i-1)) {
                covij <- cov(vars[[i]], vars[[j]])
                vcvm[i,j] <- vcvm[j,i] <- covij
            }
        }
        means <- sapply(vars, mean)
        # Generate more than needed so that we can reject those that exceed
        # the observed bounds:
        n_sims <- .n_plants + .n_plants %/% 2L
        z <- mvrnorm(n_sims, mu = means, Sigma = vcvm)
        # choose `.n_plants` rows where none exceed bounds,
        # and re-generate numbers if not enough
        good_rows_fun <- function(z, vars) {
            lgl <- z[,"mean"] >= min(vars[["mean"]]) &
                z[,"mean"] <= max(vars[["mean"]]) &
                z[,"sd"] >= min(vars[["sd"]]) & z[,"sd"] <= max(vars[["sd"]]) &
                z[,"sum"] >= min(vars[["sum"]]) & z[,"sum"] <= max(vars[["sum"]])
            return(which(lgl))
        }
        good_rows <- good_rows_fun(z, vars)
        iters <- 0L
        while (length(good_rows) < .n_plants) {
            n_sims <- n_sims + n_sims %/% 2L
            z <- mvrnorm(n_sims, mu = means, Sigma = vcvm)
            good_rows <- good_rows_fun(z, vars)
            iters <- iters + 1L
            if (iters > 5L) stop("Too many iterations to mvrnorm call.")
        }
        # choose rows that don't exceed bounds and transform back to
        # original scale:
        z <- exp(z[sample(good_rows, .n_plants),])
        args <- list(R_hat = z[,"sum"], par2 = z[,"sd"],
                     par1 = z[,"mean"] - flower_start,
                     distr_types = rep("N", .n_plants))
    } else {
        plt_idx <- sample.int(nrow(phen_fits), .n_plants,
                              replace = FALSE)
        args <- list(R_hat = phen_fits[["sum"]][plt_idx],
                     par2 = phen_fits[["sd"]][plt_idx],
                     par1 = phen_fits[["mean"]][plt_idx] - flower_start,
                     distr_types = rep("N", .n_plants))
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

    stopifnot(!is.null(nrow(other_args$z)) && nrow(other_args$z) >= 2)
    .n_plants <- nrow(other_args$z)

    for (n in names(other_args)) {
        x <- other_args[[n]]
        if (n %in% c("u", "q", "w", "add_F", "min_F_for_P")) {
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

    arg_list <- list(Y0 = rep(0, .n_plants),
                     B0 = rep(0, .n_plants),
                     # ----------*
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
                     # -------------*
                     # others:
                     P_max = rep(.n_plants / 2, .n_plants), # maximum P density
                     W = rep(0, .n_plants), # attraction to non-focal flowers (??)
                     u = 0,   # power for B/F -> P (??)
                     q = 0, # affects the strength of F -> P (??)
                     w = 0, # how quickly effects of nearby focal flowers affect P (??)
                     min_F_for_P = 0,
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
        mutate(p = factor(as.integer(p) + 1L)) |>
        filter(t > 0)

    return(sim_df)

}


time_series <- function(run_df, main = NULL, ylim_F = NULL, ylim_P = NULL,
                        plant_facet_rows = NA, P_per_flower = FALSE,
                        show_F = TRUE) {
    dd <- run_df
    if (P_per_flower) {
        dd <- dd |>
            mutate(P = ifelse((Y + B + N) == 0, P, P / (Y + B + N)))
        ylab_P <- "Pollinators / total flowers"
    } else ylab_P <- "Pollinators"
    if (show_F) {
        dd <- dd |>
            mutate(total = Y + B + N) |>
            pivot_longer(Y:total, names_to = "type", values_to = "density") |>
            mutate(type = factor(type, levels = c("total", "Y", "B", "N", "P"),
                                 labels = c("total", "yeast", "bacteria",
                                            "non-colonized", "pollinators")))
    } else {
        dd <- dd |>
            pivot_longer(Y:P, names_to = "type", values_to = "density") |>
            mutate(type = factor(type, levels = c("Y", "B", "N", "P"),
                                 labels = c("yeast", "bacteria", "non-colonized",
                                            "pollinators")))
    }
    dd <- dd |>
        mutate(id = interaction(p, type, drop = TRUE)) |>
        filter(t %% 1 == 0)
    .n_plants <- length(unique(run_df$p))
    if (.n_plants > 2) {
        .a <- 0.25
        if (.n_plants > 50) .a <- 0.1
        dd <- dd |>
            filter(type != "non-colonized")
    } else {
        .a <- 1
    }

    p1 <- dd |>
        filter(type != "pollinators") |>
        ggplot(aes(t, density)) +
        geom_hline(yintercept = 0, linewidth = 1, linetype = "22",
                   color = "gray70") +
        geom_line(aes(color = type, group = id), linewidth = 1, alpha = .a) +
        xlab("Time (days)") +
        ggtitle(main) +
        scale_y_continuous("flower-type density", limits = ylim_F) +
        scale_color_manual(NULL, values = spp_pal, guide = "none") +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5))
    p2 <- dd |>
        filter(type == "pollinators") |>
        ggplot(aes(t, (density))) +
        geom_hline(yintercept = 0, linewidth = 1, linetype = "22",
                   color = "gray70") +
        geom_line(aes(group = id), linewidth = 1, alpha = .a,
                  color = spp_pal[["pollinators"]]) +
        xlab("Time (days)") +
        scale_y_continuous(ylab_P, limits = ylim_P)
    if (!is.na(plant_facet_rows)) {
        stopifnot(is.numeric(plant_facet_rows) && length(plant_facet_rows) == 1)
        stopifnot(plant_facet_rows %% 1 == 0)
        if (.n_plants > 10) stop("I refuse to facet by plant with > 10 plants.")
        p1 <- p1 + facet_wrap(~ p, nrow = plant_facet_rows)
        p2 <- p2 + facet_wrap(~ p, nrow = plant_facet_rows)
    }
    p1 / p2
}




#' ========================================================================
#' ========================================================================
# TWO PATCHES ----
#' ========================================================================
#' ========================================================================

plt_idx <- c(20, 26)
two_patch_phen <- list(R_hat = phen_fits[["sum"]][plt_idx],
             par2 = phen_fits[["sd"]][plt_idx],
             par1 = phen_fits[["mean"]][plt_idx] - flower_start,
             distr_types = rep("N", 2))


two_patch_runs <- crossing(.q = c(0, 0.5, 1, 2, 4),
                           .u = c(0, 0.1, 1, 10)) |>
    pmap_dfr(\(.q, .u) {
        run_df <- one_run(z = 1 - diag(2L),
                          W = 50^.q * 0.5^.u,
                          q = .q,
                          u = .u,
                          w = Inf,
                          min_F_for_P = 0,
                          R_hat = two_patch_phen$R_hat,
                          add_F = 1,
                          par2 = two_patch_phen$par2,
                          par1 = two_patch_phen$par1)
        run_df |>
            mutate(q = .q, u = .u)
    }) |>
    mutate(across(q:u, factor)) |>
    select(q, u, everything())





two_patch_plots <- two_patch_runs |>
    split(~ q + u, drop = TRUE, sep = "_") |>
    map(\(x) {
        P_ymax <- ceiling(max(two_patch_runs$P) * 100) / 100
        F_ymax <- with(two_patch_runs, max(Y + B + N)) |> ceiling()
        p <- time_series(x, plant_facet_rows = 2,
                         ylim_F = c(0, F_ymax), ylim_P = c(0, P_ymax),
                         main = sprintf("q = %s | u = %s", paste(x$q[[1]]),
                                        paste(x$u[[1]]))) &
            theme(plot.title = element_text(hjust = 0.5, size = 9),
                  axis.text = element_text(size = 7),
                  axis.title = element_text(size = 8),
                  strip.text = element_blank(),
                  strip.background = element_blank())
        return(p)
    })

do.call(wrap_plots, c(two_patch_plots, list(ncol = length(levels(two_patch_runs$q)))))


cairo_pdf("_figures/two_patches.pdf", width = 12, height = 12)
do.call(wrap_plots, c(two_patch_plots, list(ncol = length(levels(two_patch_runs$q)))))
dev.off()



two_patch_summ <- two_patch_runs |>
    group_by(q, u, p) |>
    summarize(TFm = median(Y + B + N),
              Ym = median(Y / (Y + B + N)),
              Bm = median(B / (Y + B + N)),
              Nm = median(N / (Y + B + N)),
              Pm = median(log10(P)),
              .groups = "drop") |>
    rename_with(\(x) str_remove(x, "m$")) |>
    select(q, u, p, Y:B, P) |>
    # mutate(P = (P - min(P)) / diff(range(P)) * max(pmax(Y, B))) |>
    pivot_longer(Y:P, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B", "P"),
                         labels = c("yeast", "bacteria", "pollinators")))

two_patch_summ_limit_df <- two_patch_summ |>
    distinct(type, u, q) |>
    mutate(density = map(type,
                         \(tt) {
                             lgl <- two_patch_summ$type == tt
                             range(two_patch_summ$density[lgl])
                         })) |>
    unnest(density)

# two_patch_summ |>
#     ggplot(aes(as.numeric(paste(q)), density)) +
#     geom_blank(data = two_patch_summ_limit_df) +
#     geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
#     geom_point(aes(shape = p, color = type), size = 3) +
#     geom_line(aes(group = interaction(p, type), color = type, linetype = p),
#               linewidth = 0.5) +
#     facet_wrap(~ u + type, labeller = \(x) label_both(x, sep = " = "),
#                ncol = 3, scales = "free_y") +
#     scale_color_manual(NULL, values = spp_pal, guide = "none") +
#     scale_x_continuous("q")
#
# two_patch_summ |>
#     ggplot(aes(as.numeric(paste(u)), density)) +
#     geom_blank(data = two_patch_summ_limit_df) +
#     geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
#     geom_point(aes(shape = p, color = type), size = 3) +
#     geom_line(aes(group = interaction(p, type), color = type, linetype = p),
#               linewidth = 0.5) +
#     facet_wrap(~ q + type, labeller = \(x) label_both(x, sep = " = "),
#                ncol = 3, scales = "free_y") +
#     scale_color_manual(NULL, values = spp_pal, guide = "none") +
#     scale_x_continuous("u")



#' ========================================================================
#' ========================================================================
# RANDOM LANDSCAPE ----
#' ========================================================================
#' ========================================================================


np <- 100

set.seed(678902345)
xy_df <- tibble(x = runif(np), y = runif(np))
rnd_land_phen <- get_phenology(np)
dist_mat <- make_dist_mat(xy_df)

# xy_df |>
#     ggplot(aes(x, y)) +
#     geom_point() +
#     coord_equal(xlim = c(0, 1), ylim = c(0, 1))

# Example phenologies:
phen_p <- pmap_dfr(c(rnd_land_phen[-which(names(rnd_land_phen) == "distr_types")],
           list(p = 1:np)),
         \(R_hat, par2, par1, p) {
             tibble(plant = factor(p, levels = 1:np),
                    time = 0:(flower_stop - flower_start + 20),
                    r = R_hat * dnorm(time, mean = par1, sd = par2),
                    hlt = ifelse(p == np-1L, "yes", "no"))
         }) |>
    ggplot(aes(time, r, group = plant)) +
    ylab("flowering rate") +
    # geom_line(alpha = 0.25, linewidth = 1) +
    geom_line(aes(alpha = hlt, color = hlt), linewidth = 1) +
    scale_alpha_manual(values = c(0.25, 1), guide = "none") +
    scale_color_manual(values = c("black", "magenta"), guide = "none") +
    # theme(axis.title = element_blank()) +
    NULL
phen_p
# ggsave("_figures/phen_plots.png", width = 3, height = 2)

# showing effects of q:
crossing(x = 1:100, h = c(0.5, 1:3)) |>
    mutate(sum_h = map_dbl(h, \(.h) sum((1:100)^.h)), z = x^h / sum_h) |>
    ggplot(aes(x, z, color = factor(h))) +
    geom_line(linewidth = 1) +
    # theme(axis.text = element_blank(), axis.title = element_blank(),
    #       axis.ticks = element_blank()) +
    scale_color_brewer(palette = "Dark2", guide = "none")


# Takes ~ 5 sec
rnd_land_runs <- crossing(q = c(0, 1, 2),
                          u = c(0, 10, 100),
                          w = c(Inf, 100, 10, 0)) |>
    split(~ q + u + w) |>
    mclapply(\(x) {
        .q <- x[["q"]]
        .u <- x[["u"]]
        .w <- x[["w"]]
        run_df <- one_run(z = dist_mat,
                          W = 50^.q * 0.5^.u,
                          q = .q,
                          u = .u,
                          w = .w,
                          min_F_for_P = 0,
                          R_hat = rnd_land_phen$R_hat,
                          add_F = 1,
                          par2 = rnd_land_phen$par2,
                          par1 = rnd_land_phen$par1)
        run_df |>
            mutate(q = .q, u = .u, w = .w)
    }) |>
    do.call(what = bind_rows) |>
    mutate(across(q:w, factor)) |>
    select(q, u, w, everything())


# Example of one time series:
one_ts_p <- rnd_land_runs |>
    filter(q == 1, u == 10, w == 0, p == np-1L) |>
    mutate(P = P * 60) |>
    pivot_longer(Y:P, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("Y", "B", "N", "P"),
                         labels = c("yeast", "bacteria",
                                    "non-colonized", "pollinators")),
           is_P = factor(type == "pollinators")) |>
    filter(t %% 1 == 0) |>
    ggplot(aes(t, density)) +
    geom_hline(yintercept = 0, linewidth = 1, linetype = "22",
               color = "gray70") +
    geom_line(aes(color = type, linetype = is_P), linewidth = 1) +
    xlab("Time (days)") +
    scale_y_continuous("flower-type density",
                       sec.axis = sec_axis(~ .x / 60,
                                           "pollinator density")) +
    scale_color_manual(NULL, values = spp_pal, guide = "none") +
    scale_linetype_manual(NULL, values = c(1, 2), guide = "none") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
one_ts_p
# ggsave("_figures/one_ts_plot.png", width = 4, height = 2)


rnd_land_plots_sep <- rnd_land_runs |>
    split(~ q + u + w, sep = "_") |>
    map(\(x) {
        P_ymax <- ceiling(max(rnd_land_runs$P) * 100) / 100
        F_ymax <- with(rnd_land_runs, max(pmax(Y, B, N))) |> ceiling()
        p <- time_series(x, ylim_F = c(0, F_ymax), ylim_P = c(0, P_ymax),
                         show_F = FALSE,
                         main = sprintf("q = %s | u = %s | w = %s",
                                        paste(x$q[[1]]),
                                        paste(x$u[[1]]),
                                        paste(x$w[[1]]))) &
            theme(plot.title = element_text(hjust = 0.5, size = 9),
                  axis.text = element_text(size = 7),
                  axis.title = element_text(size = 8),
                  strip.text = element_blank(),
                  strip.background = element_blank())
        return(p)
    })



rnd_land_plots <- levels(rnd_land_runs$w) |>
    set_names() |>
    map(\(w) {
        idx <- which(grepl(sprintf("_%s$", w), names(rnd_land_plots_sep)))
        do.call(wrap_plots, c(rnd_land_plots_sep[idx],
                              list(ncol = length(levels(rnd_land_runs$q)))))
    })


# rnd_land_plots[["0"]]
# rnd_land_plots[["0.05"]]
# rnd_land_plots[["0.1"]]
# rnd_land_plots[["0.2"]]
# rnd_land_plots[["Inf"]]


for (n in names(rnd_land_plots)) {
    if (n == "Inf") {
        fn <- sprintf("_figures/rnd_land_ts_w=%s.pdf", n)
    } else fn <- sprintf("_figures/rnd_land_ts_w=%03i.pdf", as.integer(n))
    cairo_pdf(fn, width = 3 * length(levels(rnd_land_runs$q)),
              height = 2.5 * length(levels(rnd_land_runs$u)))
    plot(rnd_land_plots[[n]])
    dev.off()
}



rnd_land_runs |>
    mutate(f = Y + B + N) |>
    filter(f > 0) |>
    group_by(w, q, u, p) |>
    summarize(P = mean(log10(P)), .groups = "drop") |>
    getElement("P") |>
    hist()


rnd_land_microbe_summs <- rnd_land_runs |>
    filter((Y + B + N) > 0) |>
    group_by(w, q, u, p) |>
    summarize(TFm = median(Y + B + N),
              Ym = median(Y / (Y + B + N)),
              Bm = median(B / (Y + B + N)),
              Nm = median(N / (Y + B + N)),
              Pm = median(log10(P)),
              .groups = "drop") |>
    rename_with(\(x) str_remove(x, "m$")) |>
    mutate(P = ifelse(P < -6, -6, P))





rnd_land_heats <- rnd_land_microbe_summs |>
    split(~ w) |>
    map(\(x) {
        x |>
            ggplot(aes(Y, B)) +
            ggtitle(sprintf("w = %s", x$w[[1]])) +
            geom_hline(yintercept = 0, linetype = 2, color = "gray60") +
            geom_vline(xintercept = 0, linetype = 2, color = "gray60") +
            stat_bin_2d(aes(fill = after_stat(density)), binwidth = rep(0.05, 2)) +
            scale_x_continuous("Median yeast proportion", limits = c(NA, 1)) +
            scale_y_continuous("Median bacteria proportion", limits = c(NA, 1)) +
            scale_fill_viridis_c(option = "rocket", end = 0.9) +
            facet_grid(u ~ q, labeller = \(x) label_both(x, sep = " = ")) +
            theme(plot.title = element_text(hjust = 0.5)) +
            coord_equal()
    })

rnd_land_heats[["0"]]
rnd_land_heats[["0.05"]]
rnd_land_heats[["0.1"]]
rnd_land_heats[["0.2"]]
rnd_land_heats[["Inf"]]




for (n in names(rnd_land_heats)) {
    if (n == "Inf") {
        fn <- sprintf("_figures/rnd_land_heats_w=%s.pdf", n)
    } else fn <- sprintf("_figures/rnd_land_heats_w=%03i.pdf", as.integer(n))
    cairo_pdf(fn, width = 2.1 * length(levels(rnd_land_runs$q)),
              height = 2 * length(levels(rnd_land_runs$u)))
    plot(rnd_land_heats[[n]])
    dev.off()
}

#' Overall, these plots show that both types of feedbacks have an overall
#' negative effect on yeast in comparison to bacteria.
#' The exception to this is when w=0 and q=2, u has a positive effect on yeast.
#' Not sure why...


#' ... especially because mean pollinator abundance goes down with all
#' feedbacks, even when w = 0 and q = 2.




rnd_land_P_plots <- rnd_land_microbe_summs |>
    split(~ w) |>
    map(\(x) {
        x |>
            ggplot(aes(P)) +
            ggtitle(sprintf("w = %s", x$w[[1]])) +
            geom_blank(data = tibble(P = range(rnd_land_microbe_summs$P))) +
            geom_histogram(binwidth = 0.3) +
            scale_x_continuous("Median pollinators",
                               breaks = -3:0 * 2,
                               labels = parse(text = c("{} <= 10^6",
                                                       paste0("10^", -2:0 * 2)))) +
            facet_grid(u ~ q, labeller = \(x) label_both(x, sep = " = "))
    })
# rnd_land_P_plots[["Inf"]]
# rnd_land_P_plots[["100"]]
# rnd_land_P_plots[["10"]]
# rnd_land_P_plots[["0"]]


for (n in names(rnd_land_P_plots)) {
    if (n == "Inf") {
        fn <- sprintf("_figures/rnd_land_Pplots_w=%s.pdf", n)
    } else fn <- sprintf("_figures/rnd_land_Pplots_w=%03i.pdf", as.integer(n))
    cairo_pdf(fn, width = 2.5 * length(levels(rnd_land_runs$q)),
              height = 2 * length(levels(rnd_land_runs$u)))
    plot(rnd_land_P_plots[[n]])
    dev.off()
}; rm(n, fn)










#' ========================================================================
#' ========================================================================
# CLUSTERED LANDSCAPE ----
#' ========================================================================
#' ========================================================================



cls_land_phen <- rnd_land_phen

cls_radii <- c(0.01, 0.02, 0.05, 1)
w_vals <- c(2 * 0:5, Inf)



set.seed(13456978)
cls_xy_sims <- map(cls_radii,  \(.r) {
    z <- rMatClust(kappa = 25, scale = .r, mu = 10, nsim = 100)
    stopifnot(z$n >= np)
    # force all to have exactly `np` plants:
    for (j in 1:length(z)) {
        idx <- sample.int(z[[j]][["n"]], np, replace = FALSE)
        z[[j]][["x"]] <- z[[j]][["x"]][idx]
        z[[j]][["y"]] <- z[[j]][["y"]][idx]
        z[[j]][["n"]] <- np
    }
    return(z)
})
names(cls_xy_sims) <- cls_radii



cls_xy_sims |>
    imap(\(d, n) {
        tibble(x = d[[1]][["x"]], y = d[[1]][["y"]]) |>
            ggplot(aes(x, y)) +
            # ggtitle(sprintf("r = %s", n)) +
            geom_point(alpha = 0.25) +
            geom_point(shape = 1) +
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank()) +
            coord_equal(xlim = c(0, 1), ylim = c(0, 1))
    }) |>
    rev() |>
    `c`(list(nrow = 2)) |>
    do.call(what = wrap_plots)


# Takes ~3 min
t0 <- Sys.time()
cls_land_runs <- crossing(r = names(cls_xy_sims),
                          w = w_vals,
                          i = 1:100) |>
    mutate(id = 1:n()) |>
    split(~ id) |>
    mclapply(\(dd) {
        .r <- dd[["r"]]
        .w <- dd[["w"]]
        .i <- dd[["i"]]
        # .r = "0.01"; .w = 10; .i = 1
        # rm(.r, cs, dist_mat__, .q, .u, vc_mat, .W, .g_b0, mu, stdev)
        cs <- cls_xy_sims[[.r]][[.i]]
        dist_mat__ <- make_dist_mat(cs$x, cs$y)

        .q <- 2
        .u <- 1

        wW <- FALSE

        if (wW) {
            mu <- 50^.q * 0.5^.u
            stdev <- mu / 20
            vc_mat <- make_vcv_mat(dist_mat__, stdev, q = .w)
            if (.w == 0) {
                .W <- rep(mu, np)
            } else .W <- mvrnorm(1, rep(mu, np), vc_mat)
            stopifnot(all(.W > 0))
            ## hist(.W); abline(v = mu, col = "red")
            .g_b0 <- rep(0.02, np)
        } else {
            # actual mu is 0.02, but changing it to this so that vc_mat
            # is positive definite:
            mu <- 200
            stdev <- mu / 20
            vc_mat <- make_vcv_mat(dist_mat__, stdev, q = .w)
            if (.w == 0) {
                .g_b0 <- rep(mu, np)
            } else .g_b0 <- mvrnorm(1, rep(mu, np), vc_mat)
            .g_b0 <- .g_b0 * (0.02 / 200)
            stopifnot(all(.g_b0 > 0))
            ## hist(.g_b0); abline(v = mu, col = "red")
            .W <- rep(50^.q * 0.5^.u, np)
        }
        # .W <- 50^.q * 0.5^.u
        # .g_b0 <- 0.02


        run_df <- one_run(z = dist_mat__,
                          W = .W,
                          q = .q,
                          u = .u,
                          w = .w,
                          g_b0 = .g_b0,
                          min_F_for_P = 0,
                          R_hat = cls_land_phen$R_hat,
                          par2 = cls_land_phen$par2,
                          par1 = cls_land_phen$par1) |>
            group_by(p) |>
            summarize(TFm = median(Y + B + N),
                      Ym = median(Y / (Y + B + N)),
                      Bm = median(B / (Y + B + N)),
                      Nm = median(N / (Y + B + N)),
                      Pm = median(log10(P)),
                      .groups = "drop") |>
            rename_with(\(x) str_remove(x, "m$")) |>
            mutate(P = ifelse(P < -6, -6, P)) |>
            mutate(r = factor(.r, levels = names(cls_xy_sims)),
                   w = factor(.w, levels = w_vals),
                   i = factor(.i, levels = 1:100)) |>
            select(r, w, i, everything())
        return(run_df)
    }) |>
    do.call(what = bind_rows)

t1 <- Sys.time()
t1 - t0





cls_land_runs |>
    ggplot(aes(Y, B)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray60") +
    geom_vline(xintercept = 0, linetype = 1, color = "gray60") +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "gray60") +
    stat_bin_2d(aes(fill = after_stat(density)), binwidth = rep(0.05, 2)) +
    scale_x_continuous("Yeast proportion", limits = c(NA, 1)) +
    scale_y_continuous("Bacteria proportion", limits = c(NA, 1)) +
    scale_fill_viridis_c(option = "rocket", end = 0.9) +
    facet_grid(r ~ w, labeller = \(x) label_both(x, sep = " = ")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_equal()




# Average distance between any two points inside a square of side length s:
avg_pw_dist <- function(s) s * ( (2 + sqrt(2) + 5 * log(1 + sqrt(2))) / 15 )

# Simulating random points in a square of length 1:
# Takes ~2 min w 6 threads
t0 <- Sys.time()
dz <- mclapply(1:1e6L, \(i) {
    dm <- make_dist_mat(x = runif(100), y = runif(100))
    return(mean(dm[lower.tri(dm)]))
}) |>
    do.call(what = c)
t1 <- Sys.time()
t1 - t0

{
    hist(dz)
    abline(v = mean(dz), col = "red")
    abline(v = avg_pw_dist(1), col = "blue", lty = 2)
    print(mean(dz))
    print(avg_pw_dist(1))
}








#'
#'
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#'
#'
#'
#' Y_hist <- function(runs_df, xlim = NULL) {
#'     runs_df |>
#'         filter((Y + B + N) > 0) |>
#'         group_by(q, p) |>
#'         summarize(Y = mean(100 * Y / (Y + B + N)), .groups = "drop") |>
#'         ggplot(aes(Y)) +
#'         geom_histogram(bins = 20, fill = spp_pal[["yeast"]]) +
#'         # scale_x_continuous("mean % yeast") +
#'         # coord_cartesian(xlim = xlim_Y)
#'         scale_x_continuous("mean % yeast", limits = xlim) +
#'         facet_wrap(~ q, ncol = 1, scales = "free_y",
#'                    labeller = \(x) label_both(x, sep = " = "))
#' }
#'
#' P_hist <- function(run_df, xlim = NULL) {
#'     run_df |>
#'         group_by(q, p) |>
#'         summarize(P = mean(log1p(P)), .groups = "drop") |>
#'         ggplot(aes(P)) +
#'         geom_histogram(bins = 20, fill = spp_pal[["pollinators"]]) +
#'         scale_x_continuous("mean log(1 + pollinators)", limits = xlim) +
#'         facet_wrap(~ q, ncol = 1, scales = "free_y",
#'                    labeller = \(x) label_both(x, sep = " = "))
#' }
#'
#'
#'
#' animated_landscape <- function(run_df, xy_df, filename = NULL) {
#'
#'     bounds <- list(x = range(xy_df$x) + c(-1,1) * 0.1 * diff(range(xy_df$x)),
#'                    y = range(xy_df$y) + c(-1,1) * 0.1 * diff(range(xy_df$y)))
#'
#'     p <- run_df |>
#'         filter(t %% 1 == 0) |>
#'         mutate(t = as.integer(t),
#'                x = map_dbl(p, \(pp) xy_df$x[[as.integer(paste(pp))]]),
#'                y = map_dbl(p, \(pp) xy_df$y[[as.integer(paste(pp))]]),
#'                nf = Y + B + N,
#'                Y = 100 * Y / nf) |>
#'         ggplot(aes(x, y)) +
#'         geom_point(aes(size = nf, color = Y)) +
#'         scale_color_viridis_c("% yeast", limits = c(0, 100)) +
#'         scale_size("# flowers", limits = c(0, NA)) +
#'         coord_equal(xlim = bounds$x, ylim = bounds$y)
#'
#'     ap <- p + transition_time(t) +
#'         labs(title = "Day: {frame_time}")
#'
#'     if (!is.null(filename)) {
#'         anim_save(sprintf("_figures/%s.gif", filename), ap,
#'                   height = 4, width = 6, units = "in", res = 150)
#'         invisible(NULL)
#'     }
#'
#'     return(ap)
#' }
#'
#'
#'
#' #' ========================================================================
#' #' ========================================================================
#' # RANDOM LANDSCAPE ----
#' #' ========================================================================
#' #' ========================================================================
#'
#' np <- 100
#'
#' set.seed(678902345)
#' xy_df <- tibble(x = runif(np), y = runif(np))
#' dist_mat <- make_dist_mat(xy_df)
#' Y0 <- runif(100, 0, 0.5)
#' B0 <- 0.5 - Y0
#' phen <- get_phenology(np)
#'
#' # xy_df |>
#' #     ggplot(aes(x, y)) +
#' #     geom_point() +
#' #     theme(axis.title = element_blank(),
#' #           axis.text = element_blank(),
#' #           axis.ticks = element_blank()) +
#' #     coord_equal(xlim = c(0, 1), ylim = c(0, 1))
#'
#' pmap_dfr(c(phen, list(p = 1:np)), \(R_hat, sigma, mu, p) {
#'     tibble(plant = factor(p, levels = 1:np),
#'            time = 0:(flower_stop - flower_start + 20),
#'            r = R_hat * dnorm(time, mu, sigma))
#' }) |>
#'     ggplot(aes(time, r, group = plant)) +
#'     geom_line(alpha = 0.25, linewidth = 1) +
#'     theme(axis.title = element_blank())
#'
#' # showing effects of q:
#' crossing(x = 1:100, h = c(0.5, 1:3)) |>
#'     mutate(sum_h = map_dbl(h, \(.h) sum((1:100)^.h)), z = x^h / sum_h) |>
#'     ggplot(aes(x, z, color = factor(h))) +
#'     geom_line(linewidth = 1) +
#'     theme(axis.text = element_blank(), axis.title = element_blank(),
#'           axis.ticks = element_blank()) +
#'     scale_color_brewer(palette = "Dark2", guide = "none")
#'
#'
#' runs_df <- map_dfr(c(0.5, 1, 2, 3),
#'                    \(.q) {
#'                        one_run(Y0 = Y0 * 0, B0 = B0 * 0,
#'                                z = dist_mat,
#'                                W = 25^(1/.q),
#'                                q = .q,
#'                                w = 100,
#'                                R_hat = phen$R_hat,
#'                                sigma = phen$sigma,
#'                                add_F = 1,
#'                                mu = phen$mu) |>
#'                            mutate(q = factor(.q, levels = c(0.5, 1, 2, 3))) |>
#'                            dplyr::select(q, everything())
#'                    })
#' ts_ps <- runs_df |>
#'     split(~ q) |>
#'     map(\(x) {
#'         p <- time_series(x, ylim_F = c(0, 250), ylim_P = c(0, 4),
#'                          main = sprintf("q = %s", paste(x$q[[1]])))
#'         return(p)
#'     })
#'
#' do.call(wrap_plots, c(ts_ps, list(ncol = 2)))
#'
#'
#'
#' Y_p <- runs_df |> Y_hist()
#' P_p <- runs_df |> P_hist()
#'
#'
#'
#' #' ========================================================================
#' #' ========================================================================
#' # CLUSTERED LANDSCAPE ----
#' #' ========================================================================
#' #' ========================================================================
#'
#'
#'
#' cluster_radii <- c(0.02, 0.1, 1)
#' w_vals <- c(1, 10, 100)
#'
#' set.seed(13456978)
#' clust_sims <- map(cluster_radii,  \(.s) rMatClust(kappa = 25, scale = .s, mu = 10, nsim = 100))
#' for (i in 1:length(clust_sims)) stopifnot(all(map_int(clust_sims[[i]], \(x) x$n) >= np))
#' # force all to have exactly `np` plants:
#' for (i in 1:length(clust_sims)) {
#'     for (j in 1:length(clust_sims[[i]])) {
#'         idx <- sample.int(clust_sims[[i]][[j]][["n"]], np, replace = FALSE)
#'         clust_sims[[i]][[j]][["x"]] <- clust_sims[[i]][[j]][["x"]][idx]
#'         clust_sims[[i]][[j]][["y"]] <- clust_sims[[i]][[j]][["y"]][idx]
#'         clust_sims[[i]][[j]][["n"]] <- np
#'     }
#' }
#' for (i in 1:length(clust_sims)) stopifnot(all(map_int(clust_sims[[i]], \(x) x$n) == np))
#' names(clust_sims) <- cluster_radii
#'
#'
#' cl_runs_df <- map2_dfr(rep(names(clust_sims), length(w_vals)),
#'                        rep(w_vals, each = length(cluster_radii)),
#'                        \(.s, .w) {
#'                            cs <- clust_sims[[.s]][[1]]
#'                            dist_mat__ <- make_dist_mat(cs$x, cs$y)
#'                            one_run(Y0 = Y0 * 0, B0 = B0 * 0,
#'                                    z = dist_mat__,
#'                                    W = 25,
#'                                    q = 1,
#'                                    w = .w,
#'                                    R_hat = phen$R_hat,
#'                                    sigma = phen$sigma,
#'                                    mu = phen$mu) |>
#'                                mutate(s = factor(.s, levels = names(clust_sims)),
#'                                       w = factor(.w, levels = w_vals)) |>
#'                                dplyr::select(s, w, everything())
#'                        })
#'
#'
#'
#' cl_ts_ps <- cl_runs_df |>
#'     filter(s == 1) |>
#'     split(~ w) |>
#'     map(\(x) time_series(x, ylim_F = c(0, 200), ylim_P = c(0, 12),
#'                          main = sprintf("w = %s", paste(x$w[[1]]))))
#'
#'
#' do.call(wrap_plots, c(cl_ts_ps, list(ncol = 1)))
#'
#' cl_runs_df |> Y_hist2()
#' cl_runs_df |> P_hist2()
#'
#'
#' Y_hist2 <- function(cl_runs_df, xlim = NULL) {
#'     cl_runs_df |>
#'         filter((Y + B + N) > 0, s == 1) |>
#'         group_by(s, w, p) |>
#'         summarize(Y = mean(100 * Y / (Y + B + N)), .groups = "drop") |>
#'         ggplot(aes(Y)) +
#'         geom_histogram(bins = 20, fill = spp_pal[["yeast"]]) +
#'         scale_x_continuous("mean % yeast", limits = xlim) +
#'         facet_wrap(~ w, ncol = 1, labeller = \(x) label_both(x, sep = " = "))
#' }
#'
#' P_hist2 <- function(cl_runs_df, xlim = NULL) {
#'     cl_runs_df |>
#'         filter(s == 1) |>
#'         group_by(s, w, p) |>
#'         summarize(P = mean(P), .groups = "drop") |>
#'         ggplot(aes(P)) +
#'         geom_histogram(bins = 20, fill = spp_pal[["pollinators"]]) +
#'         scale_x_continuous("mean pollinator density", limits = xlim) +
#'         facet_wrap(~ w, ncol = 1, labeller = \(x) label_both(x, sep = " = "))
#' }
#'
