
# Yeq = Y[dy, L, m] = dy L Y (1 - Y - B) - m Y;
# Beq = B[db, db0, L, m] = (db0 + db L) B (1 - Y - B) - m B;
# Solve[Yeq == 0 && Beq == 0, {Y, B}]

# Y1eq = Y1'[t] = dy (Exp[-q B1] / (X + Exp[-q B1] + Exp[-q B2])) / (L0 + (Exp[-q B1] / (X + Exp[-q B1] + Exp[-q B2]))) Y1 (1 - Y1 - B1) - m Y1;
# Y2eq = Y2'[t] = dy (Exp[-q B2] / (X + Exp[-q B1] + Exp[-q B2])) / (L0 + (Exp[-q B2] / (X + Exp[-q B1] + Exp[-q B2]))) Y2 (1 - Y2 - B2) - m Y2;
# B1eq = B1'[t] = (db0 + db (Exp[-q B1] / (X + Exp[-q B1] + Exp[-q B2])) / (L0 + (Exp[-q B1] / (X + Exp[-q B1] + Exp[-q B2])))) B1 (1 - Y1 - B1) - m B1;
# B2eq = B2'[t] = (db0 + db (Exp[-q B2] / (X + Exp[-q B1] + Exp[-q B2])) / (L0 + (Exp[-q B2] / (X + Exp[-q B1] + Exp[-q B2])))) B2 (1 - Y2 - B2) - m B2;
# Solve[Y1eq == 0 && B1eq == 0 && Y2eq == 0 && B2eq == 0, {Y1, B1, Y2, B2}]






library(sweetsoursong)
library(tidyverse)





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
    d_yp = 1.1, # 1.5,  # pollinator-dependent for yeast (**)
    d_b0 = 0.3,  # pollinator-independent for bacteria (**)
    d_bp = 0.4,  # pollinator-dependent for bacteria (**)
    # -------------
    # rates of immigration from non-focal-plant sources:
    g_yp = 0, # 0.005, # pollinator-dependent for yeast (**)
    g_b0 = 0, # 0.02, # pollinator-independent for bacteria (**)
    g_bp = 0, # 0.002,  # pollinator-dependent for bacteria
    # -------------
    # pollinators:
    L_0 = 0.5,          # for P -> dispersal (**)
    q = 10,              # strength of B -> P (??)
    X = 0,              # attraction for nearby flowers other than focal plant (??)
    #
    # (**)  = value from Song et al. (submitted)
    # (??) = value should be varied bc I have no clue what to use
    #
    max_t = 1000
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
    geom_hline(yintercept = 0, linewidth = 1, linetype = "22", color = "gray70") +
    geom_line(aes(color = type, linetype = type), linewidth = 1) +
    facet_grid(run ~ p) +
    xlab("Time (days)") +
    scale_color_viridis_d(NULL, begin = 0.1, option = "H") +
    scale_linetype_manual(NULL, values = c(1, 1, 2))

