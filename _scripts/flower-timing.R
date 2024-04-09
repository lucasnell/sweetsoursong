#'
#' Quantify flowering phenology for sticky monkeyflower in Jasper Ridge
#' Biological Preserve.
#'

library(MASS)  # make sure this is loaded before tidyverse
library(moments)
library(tidyverse)
library(sweetsoursong)
library(ggtext)



phen_df <- paste0("~/Stanford_Drive/_sweetandsour/BIO47-data-2024-03-11/",
       "flower-phenology-2009.csv") |>
    read_csv(col_types = cols()) |>
    pivot_longer(-plant, names_to = "date", values_to = "n_flowers") |>
    mutate(date = as.Date(date, "%m/%d/%y"),
           plant = factor(as.integer(plant)),
           n_flowers = as.integer(round(n_flowers)),
           doy = yday(date))


phen_df |>
    ggplot(aes(date, n_flowers)) +
    geom_line(aes(color = plant)) +
    geom_point(aes(color = plant)) +
    scale_color_viridis_d(option = "A", guide = "none") +
    ylab("Number of flowers") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black"))



#' Using methods from "Skewness in bee and flower phenological distributions"
#' by Stemkovski et al. (2022)
#' (doi: 10.1002/ecy.3890)


phen_test_fits <- phen_df |>
    split(~ plant) |>
    map_dfr(\(x) {
        phen_vec <- rep(x$doy, x$n_flowers)
        skw <- skewness(phen_vec)
        skw_tst <- agostino.test(phen_vec)
        mods <- c("lognormal", "normal", "weibull") |>
            set_names() |>
            map_dbl(\(m) {
                d <- tryCatch(fitdistr(phen_vec, m),
                              error = \(e) NA,
                              warning = \(w) NA)
                if (inherits(d, "fitdistr")) {
                    return(AIC(d))
                } else return(NA)
            })
        bm <- names(mods)[which(mods == min(mods, na.rm = TRUE))]
        tibble(plant = x$plant[[1]], skew = skw, skew_p = skw_tst$p.value,
               best_mod = bm,
               dAIC_norm = mods[["normal"]] - mods[[bm]])
    }) |>
    mutate(significant = factor(skew_p < 0.05, levels = c(FALSE, TRUE),
                                labels = c("no", "yes")))

#'
#' Exactly half are skewed (15 / 30). Of these, the vast majority (13 / 15; 87%)
#' are left skewed.
#'
#' Generally, Weibull distribution fits a left-skewed phenology better, while
#' lognormal fits the right-skewed ones better.
#'
phen_test_fits |>
    ggplot(aes(best_mod, skew)) +
    geom_hline(yintercept = 0, linetype = "22", color = "gray70") +
    geom_jitter(aes(color = significant), height = 0, width = 0.1, size = 3,
                alpha = 0.5) +
    geom_text(data = tibble(best_mod = factor("lognormal"),
                            skew = c(-0.25, 0.25),
                            lab = c("left skew", "right skew"),
                            hj = c(1, 0)),
              aes(label = lab, hjust = hj),
              nudge_x = -0.5, angle = 90, vjust = 1, fontface = "italic") +
    geom_text(data = phen_test_fits |>
                  group_by(best_mod) |>
                  summarize(skew = max(phen_test_fits$skew) * 1.1,
                            n = sprintf("italic(n) == %i", n())),
              aes(label = n), parse = TRUE) +
    scale_color_manual("Significantly<br>skewed<br>(*P* < 0.05)",
                       values = c("black", "red")) +
    xlab("Best model") +
    ylab("Skew") +
    theme(legend.title = element_markdown())




phen_fits <- phen_df |>
    split(~ plant) |>
    map_dfr(\(x) {
        p <- x$plant[[1]]
        m <- phen_test_fits$best_mod[phen_test_fits$plant == p]
        phen_vec <- rep(x$doy, x$n_flowers)
        mod <- fitdistr(phen_vec, m)
        tibble(plant = p, model = m,
               par = names(mod$estimate),
               val = mod$estimate)
    })


phen_fits |>
    split(~ model) |>
    map(\(x) {
        x |>
            pivot_wider(names_from = par, values_from = val) |>
            summarize(across(c(-plant, -model), list("mean" = mean, "sd" = sd)),
                      n = n())
    })


# $lognormal
# # A tibble: 1 × 5
#   meanlog_mean meanlog_sd sdlog_mean sdlog_sd     n
#          <dbl>      <dbl>      <dbl>    <dbl> <int>
# 1         4.95      0.116     0.0756   0.0259     4
#
# $normal
# # A tibble: 1 × 5
#   mean_mean mean_sd sd_mean sd_sd     n
#       <dbl>   <dbl>   <dbl> <dbl> <int>
# 1      145.    7.72    11.9  3.20    14
#
# $weibull
# # A tibble: 1 × 5
#   shape_mean shape_sd scale_mean scale_sd     n
#        <dbl>    <dbl>      <dbl>    <dbl> <int>
# 1       15.2     3.54       151.     5.62    12
