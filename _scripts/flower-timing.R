#'
#' Quantify flowering phenology for sticky monkeyflower in Jasper Ridge
#' Biological Preserve.
#'

library(MASS)  # make sure this is loaded before tidyverse
library(moments)
library(tidyverse)

theme_set(theme_classic())


phen_df <- paste0("~/Stanford_Drive/_sweetandsour/BIO47-data-2024-03-11/",
       "flower-phenology-2009.csv") |>
    read_csv(col_types = cols()) |>
    pivot_longer(-plant, names_to = "date", values_to = "n_flowers") |>
    mutate(date = as.Date(date, "%m/%d/%y"),
           plant = factor(as.integer(plant)),
           n_flowers = as.integer(round(n_flowers)),
           doy = yday(date))

grp_phen_df <- phen_df |>
    group_by(plant) |>
    mutate(rel_flowers  = n_flowers / max(n_flowers)) |>
    group_by(date, doy) |>
    summarize(mean_flowers = mean(n_flowers),
              total_flowers = sum(n_flowers),
              rel_flowers = mean(rel_flowers),
              .groups = "drop")

phen_df |>
    ggplot(aes(date, n_flowers)) +
    geom_line(aes(color = plant)) +
    geom_point(aes(color = plant)) +
    scale_color_viridis_d(option = "A", guide = "none") +
    ylab("Number of flowers") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black"))


grp_phen_df |>
    ggplot(aes(date, mean_flowers)) +
    geom_line() +
    geom_point() +
    ylab("Number of flowers") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black"))




#' Using methods from "Skewness in bee and flower phenological distributions"
#' by Stemkovski et al. (2022)
#' (doi: 10.1002/ecy.3890)

phen_vec <- rep(grp_phen_df$doy, round(grp_phen_df$mean_flowers))
skewness(phen_vec)
agostino.test(phen_vec)
agostino.test(phen_vec, alternative = "greater")

mods <- tibble(mod = c("cauchy", "chi-squared", "exponential", "gamma",
                       "geometric", "lognormal", "logistic", "negative binomial",
                       "normal", "Poisson", "t", "weibull"),
               aic = NA_real_)

for (m in mods$mod) {
    d <- tryCatch(fitdistr(phen_vec, m), error = \(e) NA, warning = \(w) NA)
    if (inherits(d, "fitdistr")) {
        mods$aic[mods$mod == m] <- AIC(d)
    }
}

mods |>
    filter(!is.na(aic)) |>
    mutate(d_aic = aic - min(aic)) |>
    arrange(aic)


weib_mod <- fitdistr(phen_vec, "weibull")
norm_mod <- fitdistr(phen_vec, "normal")



grp_phen_df |>
    ggplot(aes(doy, mean_flowers)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_line(data = tibble(doy = seq(min(grp_phen_df$doy), max(grp_phen_df$doy),
                                      length.out = 101),
                            mean_flowers = dweibull(doy, shape = weib_mod$estimate[["shape"]],
                                         scale = weib_mod$estimate[["scale"]])) |>
                  mutate(mean_flowers = mean_flowers / max(mean_flowers) * max(grp_phen_df$mean_flowers)),
              color = "red", linewidth = 1) +
    geom_line(data = tibble(doy = seq(min(grp_phen_df$doy), max(grp_phen_df$doy),
                                      length.out = 101),
                            mean_flowers = dnorm(doy, mean = norm_mod$estimate[["mean"]],
                                         sd = norm_mod$estimate[["sd"]])) |>
                  mutate(mean_flowers = mean_flowers / max(mean_flowers) *
                             max(grp_phen_df$mean_flowers)),
              color = "gray70", linewidth = 1) +
    geom_text(data = tibble(doy = 100, lab = c("normal", "weibull"),
                            mean_flowers = c(140, 130)),
              aes(label = lab, color = lab), hjust = 0, fontface = "bold") +
    ylab("Number of flowers") +
    xlab("Day of year") +
    scale_color_manual(values = c("gray70", "red"), guide = "none")



