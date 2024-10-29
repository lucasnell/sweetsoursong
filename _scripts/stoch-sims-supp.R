
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


big_stoch_sims_file <- "_data/big-stoch-sims.rds"
# where immigration comes from other plants on the landscape so including
# it still maintains a closed system:
big_closed_stoch_sims_file <- "_data/big-stoch-closed-sims.rds"


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
# functions, objects ----
# ===========================================================================*
# ===========================================================================*


#' These are the bounds of where coexistence starts to occur (calculated
#' in Mathematica), plus the middle of these bounds (when d_b0 = 0.3 and u = 1).
#' We then moved the values that result in exclusion a bit further from the
#' border so that the outcomes don't happen extremely slowly.
d_yp__ = c("bacteria only" = 0.886308 - 0.1,
           "coexistence" = 1.17755,
           "yeast only" = 1.46879 + 0.1)





# ===========================================================================*
# ===========================================================================*
# simulations ----
# ===========================================================================*
# ===========================================================================*



if (! file.exists(big_stoch_sims_file)) {

    # Takes ~34 min with 6 threads
    set.seed(1472844374)
    stoch_sim_df <- crossing(.u = 0:10,
                             .d_yp = d_yp__,
                             .season_surv = c(0.015, 0.03, 0.06),
                             .q = c(0, 0.25, 0.5, 0.75, 1)) |>
        # Don't do this in parallel bc plant_metacomm_stoch is already
        # doing that
        pmap_dfr(\(.u, .d_yp, .season_surv, .q) {
            plant_metacomm_stoch(np = 100L, u = .u, d_yp = .d_yp,
                                 q = .q,
                                 no_immig = TRUE,
                                 n_sigma = 100,
                                 season_surv = .season_surv,
                                 # only sample last 5 seasons:
                                 burnin = 3000 * 15 / 20,
                                 # summarize by rep:
                                 summarize = "rep") |>
                mutate(u = .u, d_yp = .d_yp, season_surv = .season_surv,
                       q = .q) |>
                select(u, d_yp, season_surv, q, everything())
        }, .progress = TRUE)

    write_rds(stoch_sim_df, big_stoch_sims_file)

} else {

    stoch_sim_df <- read_rds(big_stoch_sims_file)

}










# ===========================================================================*
# ===========================================================================*
# plots ----
# ===========================================================================*
# ===========================================================================*

add_factors <- function(d) {
    if ("d_yp" %in% colnames(d) && !is.factor(d[["d_yp"]])) {
        d <- mutate(d, d_yp = factor(d_yp, levels = unname(d_yp__),
                                     labels = names(d_yp__)))
    }
    if ("q" %in% colnames(d) && !is.factor(d[["q"]])) {
        d <- mutate(d, q = factor(q))
    }
    if ("species" %in% colnames(d) && !is.factor(d[["species"]])) {
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

add_outcome_props <- function(d, .threshold = 1e-6){

    sc <- c(colnames(d)[!grepl("Y|B", colnames(d))], "outcome")
    sf <- paste("~", paste(sc, collapse = " + ")) |> as.formula()

    d |>
        mutate(outcome = case_when(meanY > .threshold & meanB > .threshold ~
                                       "coexistence",
                                   meanY > .threshold & meanB < .threshold ~
                                       "yeast only",
                                   meanY < .threshold & meanB > .threshold ~
                                       "bacteria only",
                                   TRUE ~ "extinction")) |>
        split(sf, drop = FALSE, sep = "__") |>
        imap_dfr(\(x, n) {
            if (nrow(x) > 0) {
                out <- slice(x, 1) |>
                    select(all_of(sc))
            } else {
                out <- tibble(.rows = 1)
                for (i in 1:length(sc)) out[[sc[i]]] <- str_split(n, "__")[[1]][[i]]
                for (.c in sc) out[[.c]] <- eval(call(paste0("as.", class(x[[.c]])),
                                                      out[[.c]]))
            }
            out[["prop"]] <- nrow(x) / 100
            return(out)
        }) |>
        mutate(outcome = factor(outcome, levels = names(outcome_pal)))


}




# =====================================================*
#           outcomes ----
# =====================================================*


outcome_plotter <- function(season_surv__, q__ = NULL) {

    if (is.null(q__)) {
        p <- stoch_sim_df |>
            filter(season_surv == season_surv__) |>
            select(q, u, d_yp, minY, minB, maxY, maxB, meanY, meanB) |>
            add_outcome_props() |>
            add_factors() |>
            ggplot(aes(u, prop, color = outcome)) +
            ggtitle(sprintf("s = %s", season_surv__)) +
            facet_grid(q ~ d_yp)
    } else {
        p <- stoch_sim_df |>
            add_factors() |>
            filter(season_surv == season_surv__,
                   q == q__) |>
            select(u, d_yp, minY, minB, maxY, maxB, meanY, meanB) |>
            add_outcome_props() |>
            ggplot(aes(u, prop, color = outcome)) +
            ggtitle(sprintf("s = %s, q = %s", season_surv__, q__)) +
            facet_wrap( ~ d_yp, nrow = 1)
    }
    p <- p +
        geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
        geom_point(aes(shape = outcome)) +
        geom_line() +
        scale_y_continuous("Percent of simulations",
                           limits = c(0, 1),
                           breaks = (0:4 * 25) / 100,
                           labels = c("0%", "", "50%", "", "100%")) +
        scale_x_continuous(breaks = seq(0, 10, 2.5),
                           labels = c("0", "", "5", "", "10")) +
        scale_shape_manual(NULL, values = outcome_shapes, drop = FALSE) +
        scale_color_manual(NULL, values = outcome_pal, drop = FALSE)
    return(p)
}







abundance_plotter <- function(season_surv__,
                              q__ = NULL,
                              free_y = FALSE) {

    .scales <- "fixed"
    .ylimits <- c(0, 100)
    .ybreaks <- (0:4 * 25)
    .ylabels <- c("0", "", "50", "", "100")

    if (free_y) {
        .scales <- "free_y"
        .ylimits <- NULL
        .ybreaks <- waiver()
        .ylabels <- waiver()
    }

    if (is.null(q__)) {
        dd <- stoch_sim_df |>
            filter(season_surv == season_surv__) |>
            select(q, u, d_yp, meanY, meanB) |>
            pivot_longer(meanY:meanB, names_to = "species") |>
            mutate(species = str_sub(species, 5L, 5L)) |>
            add_factors()
        dds <- dd |>
            group_by(q, u, d_yp, species) |>
            summarize(value = mean(value), .groups = "drop")
        p <- dd |>
            ggplot(aes(u, value, color = species)) +
            ggtitle(sprintf("s = %s", season_surv__)) +
            facet_grid(q ~ d_yp, scales = .scales)
    } else {
        dd <- stoch_sim_df |>
            add_factors() |>
            filter(season_surv == season_surv__,
                   q == q__) |>
            select(u, d_yp, meanY, meanB) |>
            pivot_longer(meanY:meanB, names_to = "species") |>
            mutate(species = str_sub(species, 5L, 5L)) |>
            add_factors()
        dds <- dd |>
            group_by(u, d_yp, species) |>
            summarize(value = mean(value), .groups = "drop")
        p <- dd |>
            ggplot(aes(u, value, color = species)) +
            ggtitle(sprintf("s = %s, q = %s", season_surv__, q__)) +
            facet_grid( ~ d_yp, scales = .scales)
    }

    p <- p +
        geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
        geom_point(aes(color = species), alpha = 0.1, shape = 1) +
        geom_line(data = dds, linewidth = 1) +
        scale_y_continuous("Microbial abundance", limits = .ylimits,
                           breaks = .ybreaks, labels = .ylabels) +
        scale_x_continuous(breaks = seq(0, 10, 2.5),
                           labels = c("0", "", "5", "", "10")) +
        scale_shape_manual(NULL, values = c(0, 2), drop = FALSE) +
        scale_linetype_manual(NULL, values = c(1, 1), drop = FALSE) +
        scale_color_manual(NULL, values = spp_pal, drop = FALSE,
                           aesthetics = c("color", "fill"))

    return(p)
}



outcome_plots <- map(sort(unique(stoch_sim_df$season_surv)), outcome_plotter)
abundance_plots <- map(sort(unique(stoch_sim_df$season_surv)),
                       \(s) abundance_plotter(s) +
                           theme(plot.title = element_blank()))


do.call(wrap_plots, c(outcome_plots, abundance_plots)) +
    plot_layout(byrow = TRUE,
                guides = "collect",
                ncol = length(unique(stoch_sim_df$season_surv))) # &
#     theme(legend.position = "none",
#           axis.title = element_blank(),
#           strip.text = element_blank())


# LEFT OFF ----
#' What to include in figure (all are only for s = 0.03):
#'   1. u vs outcomes (just for randomized between seasons)
#'   2. u vs microbial abundance (just for randomized between seasons)
#'   3. Something showing how non-random between seasons makes coexistence
#'      less likely at low s and more likely at high s?
#'

outcome_plotter(median(unique(stoch_sim_df$season_surv)), 0.0) /
    abundance_plotter(median(unique(stoch_sim_df$season_surv)), "random")



z <- plant_metacomm_stoch(np = 2L, u = 5, d_yp = 0.9,
                          q = 1,
                          no_immig = TRUE,
                          n_sigma = 100,
                          season_surv = 0.05,
                          max_t = 750,
                          burnin = 0,
                          summarize = "none")

tibble(t0 = c(z$Y[z$t == 150], z$B[z$t == 150]),
       t1 = c(z$Y[z$t == 150.1], z$B[z$t == 150.1]),
       species = rep(c("Y", "B"), each = sum(z$t == 150))) |>
    add_factors() |>
    ggplot(aes(t0, t1, color = species)) +
    geom_point() +
    scale_color_manual(values = spp_pal)


z |>
    filter(rep %in% 1:3) |>
    # filter(t %% 1 == 0) |>
    select(-P) |>
    pivot_longer(Y:B, names_to = "species", values_to = "density") |>
    mutate(id = interaction(species, p, rep)) |>
    add_factors() |>
    ggplot(aes(t, density, color = species)) +
    geom_line(aes(group = id), linewidth = 0.5, alpha = 0.5) +
    facet_grid(p ~ rep) +
    scale_color_manual(values = spp_pal)




logit <- gameofclones::logit
inv_logit <- gameofclones::inv_logit

par(mfrow = c(1, 4))

# cc <- 0.01
# sigma_s <- 1
q <- 0.5
Y0 <- runif(1000)
Y0[sample.int(length(Y0), length(Y0) %/% 4)] <- 0
# chi <- logit(Y0 + cc)
# eta <- chi + rnorm(length(Y0), 0, sigma_s)
# zeta <- inv_logit(eta * sd(chi) / sd(eta)) - cc
# Yt <- zeta * sum(Y0) / sum(zeta)
eps <- runif(length(Y0))
Yt <- q * Y0 + (1 - q) * eps * sum(Y0) / sum(eps)
sum(Yt); sum(Y0)


plot(Y0, Yt); hist(Yt); hist(Y0); plot(Y0[Y0 == 0], Yt[Y0 == 0])
sum(Yt); sum(Y0)
sd(Yt); sd(Y0)




sim_rand <- function(sigma_s, cc = 0.01, .fun = sd) {
    Y0 <- runif(1000, 0, 0.9)
    Y0[sample.int(length(Y0), length(Y0) %/% 4)] <- 0
    chi <- logit(Y0 + cc)
    eta <- chi + rnorm(length(Y0), 0, sigma_s)
    # zeta <- inv_logit(eta * sd(chi) / sd(eta)) - cc
    zeta <- inv_logit(eta) - cc
    zeta[zeta < 0] <- 0
    Yt <- zeta * sum(Y0) / sum(zeta)
    return(c(.fun(Y0), .fun(Yt)))
}
sim_rand2 <- function(q, .fun = sd) {
    Y0 <- runif(1000, 0, 0.9)
    Y0[sample.int(length(Y0), length(Y0) %/% 4)] <- 0
    eps <- runif(length(Y0))
    Yt <- q * Y0 + (1 - q) * mean(Y0) * (eps / sum(eps) * length(eps))
    return(c(.fun(Y0), .fun(Yt)))
}

sim_many <- function(n, .fun1, .fun2 = sd, ...) {
    sim_m <- t(replicate(n, .fun1(..., .fun = .fun2)))
    colnames(sim_m) <- c("y0", "yt")
    as_tibble(sim_m) |>
        mutate(rep = 1:n) |>
        select(rep, everything())
}

z <- crossing(sigma_s = c(0, 0.1, 0.5, 1, 5, 10),
              cc = c(0.01, 1e-3, 1e-4)) |>
    mutate(sim = map2(sigma_s, cc, \(s, .c) sim_many(1000, sim_rand, sigma_s = s, cc = .c))) |>
    unnest(sim)

z |>
    mutate(across(c(cc, sigma_s), factor)) |>
    ggplot(aes(y0, yt)) +
    geom_hline(yintercept = min(c(z$yt, z$y0)) - 0.005, color = "gray70") +
    geom_vline(xintercept = min(c(z$yt, z$y0)) - 0.005, color = "gray70") +
    geom_abline(slope = 1, intercept = 0, color = "gray70") +
    geom_point() +
    facet_grid(cc ~ sigma_s, labeller = label_both) +
    coord_equal()

z2 <- crossing(q = 0:9 / 10) |>
    mutate(sim = map(q, \(.q) sim_many(1000, sim_rand2, q = .q))) |>
    unnest(sim)

z2 |>
    mutate(across(c(q), factor)) |>
    ggplot(aes(y0, yt)) +
    geom_hline(yintercept = min(c(z2$yt, z2$y0)) - 0.005, color = "gray70") +
    geom_vline(xintercept = min(c(z2$yt, z2$y0)) - 0.005, color = "gray70") +
    geom_abline(slope = 1, intercept = 0, color = "gray70") +
    geom_point() +
    facet_wrap( ~ q) +
    coord_equal()






n <- length(Y0)
x <- cbind(Y0[1:n], runif(n)) |> unname()
y <- foo(x, sigma_s = 0.5, season_surv = 1)

plot(x[,1], y[,1]); hist(y[,1]); hist(x[,1])
colSums(x); colSums(y)

par(mfrow = c(1,1))
plot(Y0[Y0 == 0], y[Y0 == 0,1])



Rcpp::sourceCpp(code =
"
//[[Rcpp::depends(RcppArmadillo, dqrng)]]

#include <RcppArmadillo.h>
#include <pcg_random.hpp>
#include <vector>
#include <cmath>
#include <random>

// logit and inverse logit functions
inline double logit(const double& p) {
    double x = std::log(p / (1 - p));
    return x;
}
inline double inv_logit(const double& x) {
    double p = 1 / (1 + std::exp(- x));
    return p;
}
inline arma::mat logit(const arma::mat& p) {
    arma::mat x(p.n_rows, p.n_cols);
    for (size_t i = 0; i < p.n_elem; i++) x(i) = std::log(p(i) / (1 - p(i)));
    return x;
}
inline arma::mat inv_logit(const arma::mat& x) {
    arma::mat p(x.n_rows, x.n_cols);
    for (size_t i = 0; i < x.n_elem; i++) p(i) = 1 / (1 + std::exp(- x(i)));
    return p;
}

//[[Rcpp::export]]
arma::mat foo(arma::mat x, const double& sigma_s, const double& season_surv) {

    // These are 32-bit integers cast as 64-bit for downstream compatibility
    std::vector<uint64_t> tmp_seeds = as<std::vector<uint64_t>>(Rcpp::runif(4,0,4294967296));
    // Converting to two 64-bit seeds for pcg32
    std::vector<uint64_t> seeds(2);
    seeds[0] = (tmp_seeds[0]<<32) + tmp_seeds[1];
    seeds[1] = (tmp_seeds[2]<<32) + tmp_seeds[3];

    pcg32 rng;
    rng.seed(seeds[0], seeds[1]);
    std::normal_distribution<double> dist;

    double c = 0.1;

    arma::mat chi = logit(x + c);
    arma::mat eta(arma::size(x));
    arma::mat zeta(arma::size(x));
    arma::mat y(arma::size(x));

    for (size_t i = 0 ; i < eta.n_elem; i++) {
        eta(i) = chi(i) + sigma_s * dist(rng);
    }

    arma::vec zeta_mult = arma::stddev(chi).t() /
        arma::stddev(eta).t();
    for (size_t i = 0 ; i < eta.n_rows; i++) {
        for (size_t j = 0 ; j < eta.n_cols; j++) {
            zeta(i,j) = inv_logit(eta(i,j) * zeta_mult(j) - c);
        }
    }

    arma::vec x_mult = arma::sum(x).t() / arma::sum(zeta).t();
    for (size_t i = 0 ; i < zeta.n_rows; i++) {
        for (size_t j = 0 ; j < zeta.n_cols; j++) {
            y(i,j) = season_surv * zeta(i,j) * x_mult(j);
        }
    }
    return y;
}")
