
#'
#' Stochastic simulations for supplemental figures.
#'

library(sweetsoursong)
library(tidyverse)
library(patchwork)
library(ggtext)
library(grid)
library(RcppParallel)
setThreadOptions(numThreads = max(defaultNumThreads() - 2L, 1L))

if (file.exists(".Rprofile")) source(".Rprofile")

# sims of outcomes / abundances:
main_stoch_sims_file <- "_data/big-stoch-sims.rds"
# sims of gain / loss:
gain_loss_stoch_sims_file <- "_data/big-stoch-gain-loss-sims.rds"
# sims of mutual invasibility:
mutual_inv_stoch_sims_file <- "_data/big-stoch-mutual-inv-sims.rds"


outcome_pal <- c("coexistence" = "#008B00",
                 "yeast only" = "#FFCC33",
                 "bacteria only" = "#333399",
                 "extinction" = "gray60")
outcome_shapes <- c("coexistence" = 19,
                    "yeast only" = 19,
                    "bacteria only" = 19,
                    "extinction" = 4)

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "gray60")


# ===========================================================================*
# ===========================================================================*
# data functions, d_b0 ----
# ===========================================================================*
# ===========================================================================*


#' "coexistence" is the boundary where coexistence occurs at u=0.
#' The other two are on either side of that boundary and produce competitive
#' exclusion.
d_b0__ = c("bacteria only" = 0.4,
           "coexistence" = 0.3,
           "yeast only" = 0.2)




add_factors <- function(d, .exclude = NULL) {
    if ("d_b0" %in% colnames(d) && !is.factor(d[["d_b0"]]) &&
        ! "d_b0" %in% .exclude) {
        d <- mutate(d, d_b0 = factor(d_b0, levels = unname(d_b0__),
                                     labels = names(d_b0__)))
    }
    if ("q" %in% colnames(d) && !is.factor(d[["q"]]) &&
        ! "q" %in% .exclude) {
        d <- mutate(d, q = factor(q))
    }
    if ("u" %in% colnames(d) && !is.factor(d[["u"]]) &&
        ! "u" %in% .exclude) {
        d <- mutate(d, u = factor(u))
    }
    if ("species" %in% colnames(d) && !is.factor(d[["species"]]) &&
        ! "species" %in% .exclude) {
        if (all(d$species %in% c("Y", "B"))) {
            d <- mutate(d, species = factor(species, levels = c("Y", "B"),
                                            labels = c("yeast", "bacteria")))
        } else if (all(d$species %in% c("yeast", "bacteria"))) {
            d <- mutate(d, species = factor(species,
                                            levels = c("yeast", "bacteria")))
        } else if (all(d$species %in% c("Y", "B", "P"))) {
            d <- mutate(d, species = factor(species, levels = c("Y", "B", "P"),
                                            labels = c("yeast", "bacteria",
                                                       "pollinators")))
        } else if (all(d$species %in% c("yeast", "bacteria", "pollinators"))) {
            d <- mutate(d, species = factor(species,
                                            levels = c("yeast", "bacteria",
                                                       "pollinators")))
        } else {
            stop("strange values in d$species")
        }
    }
    return(d)
}

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




# ===========================================================================*
# ===========================================================================*
# even sims ----
# ===========================================================================*
# ===========================================================================*



if (! file.exists(main_stoch_sims_file)) {

    # Takes ~10 min with 6 threads
    set.seed(1472844374)
    stoch_sim_df <- crossing(.u = 0:10,
                             .d_b0 = d_b0__,
                             .q = c(0, 0.25, 0.5, 0.75, 0.95)) |>
        # Don't do this in parallel bc plant_metacomm_stoch is already
        # doing that
        pmap_dfr(\(.u, .d_b0, .q) {
            plant_metacomm_stoch(np = 100L, u = .u, d_b0 = .d_b0,
                                 q = .q,
                                 no_immig = TRUE,
                                 n_sigma = 100,
                                 # only sample last season:
                                 burnin = 3000 - 150,
                                 # summarize by rep:
                                 summarize = "rep") |>
                mutate(u = .u, d_b0 = .d_b0, q = .q) |>
                select(u, d_b0, q, everything())
        }, .progress = TRUE)

    write_rds(stoch_sim_df, main_stoch_sims_file)

} else {

    stoch_sim_df <- read_rds(main_stoch_sims_file)

}




# =====================================================*
#           mutual inv. sims ----
# =====================================================*



# if (! file.exists(mutual_inv_stoch_sims_file)) {
#
#     # Takes ~26 min with 6 threads
#     t0 <- Sys.time()
#     set.seed(1771251461)
#     inv_sim_df <- crossing(.u = 0:10,
#                            .d_b0 = d_b0__,
#                            .q = sort(unique(stoch_sim_df$q)),
#                            .rare_spp = c("yeast", "bacteria")) |>
#         # Don't do this in parallel bc plant_metacomm_stoch is already
#         # doing that
#         pmap_dfr(\(.u, .d_b0, .q, .rare_spp) {
#             if (.rare_spp == "yeast") {
#                 .Y0 <- c(rep(0.5, 1L), rep(0,   99L))
#                 .B0 <- c(rep(0,   1L), rep(0.5, 99L))
#             } else {
#                 .B0 <- c(rep(0.5, 1L), rep(0,   99L))
#                 .Y0 <- c(rep(0,   1L), rep(0.5, 99L))
#             }
#             z <- plant_metacomm_stoch(np = 100L, u = .u, d_b0 = .d_b0, q = .q,
#                                       Y0 = .Y0, B0 = .B0,
#                                       no_immig = TRUE,
#                                       n_sigma = 100,
#                                       max_t = 3000,
#                                       # only calculate abundance based on the
#                                       # last half of the last season:
#                                       burnin = 3000 - 75 + 0.1) |>
#                 select(-P) |>
#                 pivot_longer(Y:B, names_to = "species", values_to = "density")
#             x <- z |>
#                 filter(t == max(t)) |>
#                 mutate(occup = density > 0) |>
#                 group_by(rep, species) |>
#                 summarize(occup = sum(occup), .groups = "drop")
#             y <- z |>
#                 group_by(rep, species, t) |>
#                 summarize(density = sum(density), .groups = "drop") |>
#                 group_by(rep, species) |>
#                 summarize(density = mean(log1p(density)), .groups = "drop")
#             bind_cols(x, density = y$density) |>
#                 mutate(u = .u, d_b0 = .d_b0, q = .q, rare_spp = .rare_spp) |>
#                 add_factors(.exclude = NULL) |>
#                 mutate(rare_spp = factor(rare_spp, levels = c("yeast",
#                                                               "bacteria"))) |>
#                 select(u, d_b0, q, rare_spp, everything())
#         }, .progress = TRUE)
#
#     write_rds(inv_sim_df, mutual_inv_stoch_sims_file)
#
#     t1 <- Sys.time()
#     t1 - t0; rm(t0, t1)
#
#
# } else {
#
#     inv_sim_df <- read_rds(mutual_inv_stoch_sims_file)
#
# }



# LEFT OFF ----

calc_P <- function(u, species, rare_sp) {

    # for patch i, equation is: n * (1 - B[i])^u / (2 * sum((1 - B[j])^u))
    # Below,
    # Bi = (1 - B[i])^u
    # sumBj = sum((1 - B[j])^u)
    # n = 100
    #
    # I also assume the rare species occupies one plant with density = 0.5
    # (the other with density = 0)
    # The common species occupies 99 plants of the 100 total plants

    if (length(species) == 1) species <- rep(species, length(u))
    if (length(rare_sp) == 1) rare_sp <- rep(rare_sp, length(u))

    stopifnot(length(u) == length(species))
    stopifnot(length(u) == length(rare_sp))

    Bi <- numeric(length(u))
    spY <- species == "yeast"
    Bi[spY] <- 1
    Bi[!spY] <- 0.5^(u[!spY])

    sumBj <- numeric(length(u))
    rsY <- rare_sp == "yeast"
    sumBj[rsY] <- 99 * 0.5^(u[rsY]) + 1
    sumBj[!rsY] <- 99 + 0.5^(u[!rsY])

    P <- 100 * Bi / (2 * sumBj)

    return(P)
}


calc_dNdt <- function(P, species, d_yp, d_b0,
                      d_bp = formals(plant_metacomm_stoch)[["d_bp"]],
                      L0 = 0.5, Y = 0.5, B = 0.5,
                      N = 0.5, m = 0.1) {

    if (length(species) == 1) species <- rep(species, length(P))
    if (length(d_yp) == 1) d_yp <- rep(d_yp, length(P))
    if (length(d_b0) == 1) d_b0 <- rep(d_b0, length(P))
    if (length(d_bp) == 1) d_bp <- rep(d_bp, length(P))

    stopifnot(length(P) == length(species))
    stopifnot(length(P) == length(d_yp))
    stopifnot(length(P) == length(d_b0))
    stopifnot(length(P) == length(d_bp))
    dNdt <- numeric(length(P))
    isY <- species == "yeast"
    .delta <- d_yp[isY] * (P[isY] / (L0 + P[isY])) * Y
    .gamma <- (d_b0[!isY] + d_bp[!isY] * (P[!isY] / (L0 + P[!isY]))) * B
    dNdt[isY] <- .delta * N - m * Y
    dNdt[!isY] <- .gamma * N - m * B
    return(dNdt)
}



invasion_calcs <- crossing(u = 0:10,
                           d_b0 = c(0.24, 0.3, 0.36),
                           d_yp = c(0.8, 1, 1.2),
                           species = c("yeast", "bacteria"),
                           rare_sp = c("yeast", "bacteria")) |>
    mutate(P = calc_P(u, species, rare_sp),
           dN = calc_dNdt(P, species, d_yp, d_b0))


# invasion_calcs |>
#     filter(species == rare_sp) |>
#     # These don't affect P:
#     filter(d_b0 == median(d_b0), d_yp == median(d_yp)) |>
#     ggplot(aes(u, P, color = species)) +
#     geom_point(size = 2) +
#     geom_line(linewidth = 1) +
#     scale_color_manual(NULL, values = spp_pal, guide = "none")
#
# invasion_calcs |>
#     filter(species == rare_sp) |>
#     select(-rare_sp, -dN) |>
#     # These don't affect P:
#     filter(d_b0 == median(d_b0), d_yp == median(d_yp)) |>
#     pivot_wider(names_from = species, values_from = P) |>
#     mutate(rel_P = yeast / bacteria) |>
#     ggplot(aes(u, rel_P)) +
#     geom_point(size = 2) +
#     geom_line(linewidth = 1)
#
#
# invasion_calcs |>
#     filter(species == rare_sp) |>
#     ggplot(aes(u, dN, color = species)) +
#     geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
#     geom_point(size = 2) +
#     geom_line(linewidth = 1) +
#     scale_color_manual(NULL, values = spp_pal, guide = "none") +
#     facet_grid(d_b0 ~ d_yp, labeller =
#                    labeller(d_yp = \(x) { lapply(x, \(z) paste("*d<sub>yp</sub>* =", z)) },
#                             d_b0 = \(x) { lapply(x, \(z) paste("*d<sub>b0</sub>* =", z)) })) +
#     theme(strip.text = element_markdown())
#
#
# invasion_calcs |>
#     select(-P) |>
#     pivot_wider(names_from = species, values_from = dN) |>
#     mutate(inv_rel_dN = case_when(rare_sp == "yeast" ~ yeast / bacteria,
#                               TRUE ~ bacteria / yeast)) |>
#     ggplot(aes(u, log2(inv_rel_dN), color = rare_sp)) +
#     geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
#     geom_point(size = 2) +
#     geom_line(linewidth = 1) +
#     scale_color_manual(NULL, values = spp_pal, guide = "none") +
#     facet_grid(d_b0 ~ d_yp, labeller =
#                    labeller(d_yp = \(x) { lapply(x, \(z) paste("*d<sub>yp</sub>* =", z)) },
#                             d_b0 = \(x) { lapply(x, \(z) paste("*d<sub>b0</sub>* =", z)) })) +
#     theme(strip.text = element_markdown())



# d_yp = 2 * d_b0 + 0.4



# Takes ~5 min
growth_df <- crossing(.u = c(0, 1, 10),
                      # .d_b0 = c(0.1, 0.2, 0.3, 0.4),
                      .d_b0 = 0.3,
                      # .d_yp = c(0.6, 0.8, 1, 1.2),
                      .d_bp = c(0.3, 0.4, 0.5),
                      .d_yp = 2 * .d_b0 + .d_bp,
                      .rare_sp = c("yeast", "bacteria")) |>
    pmap_dfr(\(.u, .d_b0, .d_bp, .d_yp, .rare_sp) {
        if (.rare_sp == "yeast") {
            .Y0 <- c(rep(0.5, 1L), rep(0,   99L))
            .B0 <- c(rep(0,   1L), rep(0.5, 99L))
        } else {
            .B0 <- c(rep(0.5, 1L), rep(0,   99L))
            .Y0 <- c(rep(0,   1L), rep(0.5, 99L))
        }
        plant_metacomm_stoch(np = 100L, n_reps = 100L,
                             u = .u,
                             Y0 = .Y0, B0 = .B0,
                             d_b0 = .d_b0,
                             d_bp = .d_bp,
                             d_yp = .d_yp,
                             q = 0.5,
                             no_immig = TRUE,
                             n_sigma = 100,
                             max_t = 3000,
                             # burnin = 3000 - 600,
                             # summarize = "time",
                             # save_every = 1)
                             begin_end = TRUE) |>
            filter(round(t %% 150, 1) == 0.1) |>
            mutate(dY = calc_dNdt(P, "yeast", d_yp = .d_yp, d_b0 = .d_b0,
                                  Y = Y, B = B, N = 1 - Y - B),
                   dB = calc_dNdt(P, "bacteria", d_yp = .d_yp, d_b0 = .d_b0,
                                  Y = Y, B = B, N = 1 - Y - B)) |>
            group_by(rep, t) |>
            summarize(p_inv = ifelse(.rare_sp == "yeast", mean(dY > dB), mean(dB > dY)),
                      .groups = "drop") |>
            mutate(u = .u, d_b0 = .d_b0, d_bp = .d_bp, d_yp = .d_yp, rare_sp = .rare_sp)
    }, .progress = TRUE)


b_disp <- "d_bp"

growth_ylims <- growth_df |>
    group_by(t, d_yp, .data[[b_disp]], rare_sp, u) |>
    summarize(p_inv = mean(p_inv), .groups = "drop") |>
    filter(p_inv > 0, p_inv < 1) |>
    getElement("p_inv") |> range()




growth_plots <- map(c(0, 1, 10), \(.u) {
    growth_df |>
        filter(u == .u) |>
        ggplot(aes(t, logit(p_inv), color = rare_sp)) +
        ggtitle(sprintf("*u* = %i", .u)) +
        # geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
        geom_hline(yintercept = logit(1/100), linetype = 1, color = "gray80") +
        # geom_point(shape = 1, size = 2, alpha = 0.5) +
        geom_line(data = growth_df |> filter(u == .u) |>
                      group_by(t, d_yp, .data[[b_disp]], rare_sp) |>
                      summarize(p_inv = mean(p_inv), .groups = "drop"),
                  linewidth = 1) +
        scale_color_manual(values = spp_pal, guide = "none") +
        scale_y_continuous("Proportion flowers where invader wins",
                           limits = logit(growth_ylims)) +
        facet_grid(vars(.data[[b_disp]]), vars(d_yp), labeller =
                       labeller(d_yp = \(x) { lapply(x, \(z) paste("*d<sub>yp</sub>* =", z)) },
                                d_b0 = \(x) { lapply(x, \(z) paste("*d<sub>b0</sub>* =", z)) },
                                d_bp = \(x) { lapply(x, \(z) paste("*d<sub>bp</sub>* =", z)) })) +
        theme(strip.text = element_markdown(),
              plot.title = element_markdown())
})


do.call(wrap_plots, growth_plots) +
    plot_layout(nrow = 1)

growth_df |>
    group_by(d_yp, .data[[b_disp]], rare_sp, u) |>
    summarize(p_inv = mean(p_inv), .groups = "drop") |>
    ggplot(aes(factor(u), logit(p_inv), color = rare_sp)) +
    geom_hline(yintercept = logit(1/100), linetype = 1, color = "gray80") +
    geom_point(size = 2) +
    geom_line(aes(x = as.numeric(factor(u))), linewidth = 1) +
    scale_color_manual(values = spp_pal, guide = "none") +
    scale_y_continuous("Proportion flowers where invader wins",
                       limits = logit(growth_ylims)) +
    facet_grid(vars(.data[[b_disp]]), vars(d_yp), labeller =
                   labeller(d_yp = \(x) { lapply(x, \(z) paste("*d<sub>yp</sub>* =", z)) },
                            d_b0 = \(x) { lapply(x, \(z) paste("*d<sub>b0</sub>* =", z)) },
                            d_bp = \(x) { lapply(x, \(z) paste("*d<sub>bp</sub>* =", z)) })) +
    theme(strip.text = element_markdown(),
          plot.title = element_markdown())



z |>
    filter(round(t %% 150, 1) == 0.1) |>
    mutate(dY = calc_dNdt(P, "yeast", d_yp = .d_yp, d_b0 = .d_b0,
                          Y = Y, B = B, N = 1 - Y - B),
           dB = calc_dNdt(P, "bacteria", d_yp = .d_yp, d_b0 = .d_b0,
                          Y = Y, B = B, N = 1 - Y - B)) |>
    ggplot(aes(t, dY - dB)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
    geom_point(alpha  = 0.5, shape = 1) +
    stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE) +
    facet_wrap(~ rep, nrow = 2)


z |>
    ggplot(aes(t, p_inv)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
    geom_point(size = 2) +
    # geom_line(linewidth = 1) +
    # facet_wrap(~ rep, nrow = 2) +
    ylim(0, 1)




z |>
    select(rep, t, sumY, sumB) |>
    rename(Y = sumY, B = sumB) |>
    pivot_longer(Y:B, names_to = "species", values_to = "density") |>
    add_factors() |>
    ggplot(aes(t, density, color = species)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
    geom_line(linewidth = 1) +
    ylab("Abundance") +
    facet_wrap(~ rep, nrow = 2) +
    scale_color_manual(NULL, values = spp_pal, guide = "none")

z |>
    select(rep, t, occupY, occupB) |>
    rename(Y = occupY, B = occupB) |>
    pivot_longer(Y:B, names_to = "species", values_to = "density") |>
    add_factors() |>
    ggplot(aes(t, density, color = species)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
    geom_line(linewidth = 1) +
    ylab("Occupancy") +
    facet_wrap(~ rep, nrow = 2) +
    scale_color_manual(NULL, values = spp_pal, guide = "none")









# ===========================================================================*
# ===========================================================================*
# plot functions ----
# ===========================================================================*
# ===========================================================================*



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
        .facet <- facet_fxn(vars(.data[[.f_var]]), vars(d_b0),
                             labeller = labeller(!!.f_var := f_labeller,
                                                 d_b0 = label_value))
    } else {
        if (is.null(facet_fxn)) facet_fxn <- facet_wrap
        .facet <- facet_fxn( ~ d_b0, nrow = 1)
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
        filter(.data[[f_var]] %in% f_lvls) |>
        select(q, u, d_b0, maxY, maxB) |>
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
        .title <- ""
        if (length(f_lvls) == 1) {
            .title <- paste0(.title, sprintf(", *%s* = %s", f_var, f_lvls))
        }
        p <- p + ggtitle(.title)
    }

    return(p)
}





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
                                        ncol = length(unique(stoch_sim_df$d_b0)),
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
        filter(.data[[f_var]] %in% f_lvls) |>
        select(q, u, d_b0, {{ yvar }}, {{ bvar }}) |>
        pivot_longer({{ yvar }}:{{ bvar }}, names_to = "species") |>
        mutate(species = str_sub(species, nchar(species), nchar(species))) |>
        add_factors(.exclude = x_var)
    dds <- dd |>
        group_by(q, u, d_b0, species) |>
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
        .title <- ""
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






invasion_plotter <- function(x_var, y_var, double_facet, plot_rare = TRUE) {

    y_var <- match.arg(y_var, c("occup", "density"))

    f_var <- ifelse(x_var == "u", "q", "u")

    lgl_fun <- ifelse(plot_rare, `==`, `!=`)

    if (x_var == "u") {
        x_axis <- scale_x_continuous(paste("Strength of microbe--pollinator",
                                           "effect (*u*)"),
                                     breaks = seq(0, 10, 2.5),
                                     labels = c("0", "", "5", "", "10"))
        dd <- inv_sim_df
    } else {
        x_axis <- scale_x_continuous("Between-season determinism (*q*)",
                                     breaks = as.numeric(levels(inv_sim_df$q)),
                                     labels = c("0", "", "0.5", "", "0.95"))
        dd <- inv_sim_df |> filter(u %in% c(0, 1, 4, 7, 10))
    }
    if (y_var == "occup") {
        y_axis <- scale_y_continuous("Last-season occupancy",
                                     breaks = log1p(c(0, 3^(0:3))),
                                     labels = c(0, 3^(0:3)))
        y0 <- 1
    } else {
        y_axis <- scale_y_continuous("Last-season abundance",
                                     breaks = log1p(2^(0:3 * 2L - 1L)),
                                     labels = 2^(0:3 * 2L - 1L))
        y0 <- 0.5
    }

    dd <- dd |>
        mutate(occup = log1p(occup)) |> # density is already transformed
        group_by(u, d_b0, q, rare_spp, species) |>
        summarize(lo = quantile(.data[[y_var]], 0.2),
                  hi = quantile(.data[[y_var]], 0.8),
                  !!y_var := mean(.data[[y_var]]),
                  .groups = "drop") |>
        filter(lgl_fun(rare_spp, species)) |>
        mutate(increasing = factor(.data[[y_var]] > log1p(y0))) |>
        mutate(!!x_var := as.numeric(paste(.data[[x_var]])))

    if (double_facet) {
        facet <- facet_grid(vars(.data[[f_var]]), vars(d_b0),
                            labeller = labeller(!!f_var := \(x) lapply(
                                x, \(z) paste0("*", f_var, "* = ", z)),
                                       d_b0 = label_value))
    } else {
        dd <- dd |>
            filter(.data[[f_var]] == median(as.numeric(levels(.data[[f_var]]))))
        facet <- facet_grid(~ d_b0)
    }

    dd |>
        ggplot(aes(.data[[x_var]], .data[[y_var]], color = species)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_hline(yintercept = log1p(y0), linetype = "22", color = "gray70") +
        geom_point(aes(shape = increasing), size = 2) +
        geom_line(linewidth = 1, show.legend = TRUE) +
        x_axis +
        y_axis +
        facet +
        scale_color_manual(NULL, values = spp_pal) +
        scale_shape_manual(NULL, values = c("TRUE" = 19, "FALSE" = 1),
                           guide = "none") +
        theme(axis.title.x = element_markdown(margin = margin(0,0,0,t=12)),
              strip.text = element_markdown())

}






# =====================================================*
#           creating plots ----
# =====================================================*




#' #'
#' #' Supplemental figures with all combos:
#' #'
#' for (x in c("u", "q")) {
#'     for (y in c("occup", "density")) {
#'         p <- invasion_plotter(x, y, double_facet = TRUE)
#'         save_plot(sprintf("_figures/supp-%s-%s-invasion.pdf", x, y), p, 4.5, 6)
#'     }
#' }; rm(p, x, y)


supp_u_out_abund_p <-
    (outcome_plotter("u", drop_ext = FALSE) +
         theme(plot.title = element_blank())) +
    (abundance_plotter("u") + theme(plot.title = element_blank())) +
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
supp_u_out_abund_p


# save_plot("_figures/supp-u-out-abund.pdf", supp_u_out_abund_p, 9, 6)

# LEFT OFF #2 -----
# - RE-MAKE (BC OF ACCIDENTAL DELETION) FIGURE `supp-u-out-abund.pdf`
# - MAKE SURE THIS FIGURE MAKES SENSE NOW THAT YOU'RE CHANGING `d_b0`
#   INSTEAD OF `d_yp`
# - ALSO ADJUST THE ORDER OF THESE PLOTS TO MATCH THE MANUSCRIPT


#' #'
#' #' Main text figure with `q = 0.5`:
#' #'
#' main_plots <- list(outcome_plotter("u", facet_q = 0.5, facet_u = 4),
#'                    abundance_plotter("u", facet_q = 0.5, facet_u = 4),
#'                    invasion_plotter("u", "occup", double_facet = FALSE)) |>
#'     map(\(p) {
#'         p +
#'             theme(plot.title = element_blank(),
#'                   legend.position = "none",
#'                   panel.spacing.x = unit(1, "lines"),
#'                   axis.title.y = element_blank(),
#'                   axis.title.x = element_blank(),
#'                   strip.text = element_blank())
#'     }) |>
#'     do.call(what = wrap_plots) +
#'     plot_layout(ncol = 1)
#'
#' save_plot("_figures/u-stochastic-main.pdf", main_plots, 4, 5)
