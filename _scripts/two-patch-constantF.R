

library(sweetsoursong)
library(tidyverse)
library(parallel)

options("mc.cores" = max(1L, detectCores()-2))



equil_run <- function(other_args, .precision = 1e-4) {

    stopifnot(inherits(other_args, "list"))
    stopifnot(all(map_int(other_args, length) == 1))

    args <- list(dt = 0.1,
                 max_t = 1000,
                 m = rep(0.1, 2),
                 # -------------
                 # rates of immigration from non-focal-plant sources:
                 g_yp = rep(0,2), # 0.005, # pollinator-dependent for yeast (**)
                 g_b0 = rep(0,2), # 0.02, # pollinator-independent for bacteria (**)
                 g_bp = rep(0,2), # 0.002,  # pollinator-dependent for bacteria
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
        if (n %in% c("P_total", "q", "X")) next
        if (length(other_args[[n]]) == 1) other_args[[n]] <- rep(other_args[[n]], 2)
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
                         q = c(0, 1, 5, 10),              # strength of B -> P (??)
                         max_t = 10e3L) |>
        (\(x) split(x, 1:nrow(x)))() |>
        map(as.list) |>
        mclapply(equil_run) |>
        do.call(what = bind_rows)
    write_rds(equil_df, "_data/two-patch-constF-equil.rds")
} else {
    equil_df <- read_rds("_data/two-patch-constF-equil.rds")
}




outcome_pal <- c("coexist" = "#4477AA",
                 "coexist - within patch" = "#66CCEE",
                 "yeast only" = "#CCBB44",
                 "bacteria only" = "#AA3377")

equil_df |>
    filter(L_0 == 0.5, d_bp == 0.4) |>
    mutate(Yprop = (Y1 + Y2) / (Y1 + Y2 + B1 + B2),
           q = factor(q, levels = sort(unique(q)),
                      labels = sprintf("q = %i", sort(unique(q)))),
           outcome = case_when(Yprop < 1e-3 ~ "bacteria only",
                               Yprop > (1-1e-3) ~ "yeast only",
                               Y1 >= 1e-3 & Y1 <= (1-1e-3) &
                                   Y2 >= 1e-3 & Y2 <= (1-1e-3) ~ "coexist - within patch",
                               TRUE ~ "coexist")) |>
    ggplot(aes(d_yp, d_b0, fill = outcome)) +
    geom_raster() +
    scale_fill_manual(values = outcome_pal) +
    xlab("Yeast pollinator-dependent dispersal") +
    ylab("Bacteria pollinator-independent dispersal") +
    facet_grid(run ~ q)







