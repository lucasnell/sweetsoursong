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
    group_by(plant) |>
    summarize(n_flowers = sum(n_flowers)) |>
    ggplot(aes(n_flowers)) +
    geom_histogram(bins  = 10)


phen_df |>
    ggplot(aes(date, n_flowers)) +
    geom_line(aes(color = plant)) +
    geom_point(aes(color = plant)) +
    scale_color_viridis_d(option = "A", guide = "none", end = 0.9) +
    ylab("Number of flowers") +
    facet_wrap(~ plant, nrow = 6) +
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
        doy0 <- min(x$doy[x$n_flowers > 0])
        tibble(plant = x$plant[[1]], skew = skw, skew_p = skw_tst$p.value,
               best_mod = bm,
               dAIC_norm = mods[["normal"]] - mods[[bm]],
               first_flower = doy0)
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

#' Plants that were already flowering when sampling started:
ff_plants <- phen_test_fits |>
    filter(first_flower == min(phen_df$doy)) |>
    getElement("plant")

phen_df |>
    mutate(ff = factor(plant %in% ff_plants, levels = c(TRUE, FALSE),
                       labels = c("flowering @ t=0", "no flowers @ t=0"))) |>
    group_by(plant) |>
    mutate(rel_flowers = n_flowers / max(n_flowers)) |>
    ungroup() |>
    ggplot(aes(date, rel_flowers)) +
    geom_line(aes(color = ff)) +
    geom_point(aes(color = ff)) +
    scale_color_manual(NULL, values = c("red", "black")) +
    ylab("Relative number of flowers") +
    facet_wrap(~ plant, nrow = 6) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black"))

phen_test_fits |> filter(significant == "yes") |> mutate(skew_type = ifelse(skew < 0, "left", "right"))


phen_test_fits2 <- filter(phen_test_fits, ! plant %in% ff_plants)

phen_test_fits |>
    ggplot(aes(first_flower, skew)) +
    geom_point() +
    stat_smooth(method = "lm", se = TRUE, formula = y ~ x + I(x^2))

phen_test_fits2 |>
    ggplot(aes(first_flower, skew)) +
    geom_point() +
    stat_smooth(method = "lm", se = TRUE, formula = y ~ x)


lm(skew ~ first_flower, phen_test_fits) |> summary()
lm(skew ~ first_flower + I(first_flower^2), phen_test_fits) |> summary()

lm(skew ~ first_flower, phen_test_fits2) |> summary()
lm(skew ~ first_flower + I(first_flower^2), phen_test_fits2) |> summary()




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



phen_fits |>
    ggplot(aes(val, 1)) +
    geom_jitter(shape = 1, height = 0.2, width = 0) +
    facet_wrap(~ model + par, scales = "free_x", ncol = 2) +
    scale_y_continuous(NULL, limits = c(0.5, 1.5)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())



# If we forced them all to have a normal distribution:

norm_fits <- phen_df |>
    split(~ plant) |>
    map_dfr(\(x) {
        p <- x$plant[[1]]
        phen_vec <- rep(x$doy, x$n_flowers)
        mod <- fitdistr(phen_vec, "normal")
        tibble(plant = p,
               sum = sum(x$n_flowers),
               mean = mod$estimate[["mean"]],
               sd = mod$estimate[["sd"]])
    })

write_rds(norm_fits, "_data/norm-phen-fits.rds")

norm_fits |>
    print(n = 30)

# # A tibble: 30 × 4
#    plant   sum  mean    sd
#    <fct> <int> <dbl> <dbl>
#  1 1       366  153. 12.3
#  2 2       251  155. 12.0
#  3 3      1294  155. 11.9
#  4 4        36  147.  7.41
#  5 5       124  156.  9.23
#  6 6      1032  134. 11.7
#  7 7      2272  135. 20.0
#  8 8      1617  131. 13.3
#  9 9       993  140. 13.0
# 10 10     1041  144.  8.87
# 11 11      914  145. 11.4
# 12 12      762  150.  8.58
# 13 13      417  156.  7.51
# 14 14     1575  125. 12.4
# 15 15      452  148.  9.23
# 16 16      400  152.  8.33
# 17 17     1396  154.  8.76
# 18 18      888  143. 12.9
# 19 19      468  148. 11.0
# 20 20      552  152. 10.3
# 21 21      473  143. 14.3
# 22 22     1221  155. 10.5
# 23 23      433  139. 13.3
# 24 24     1014  138. 14.2
# 25 25      324  146. 12.6
# 26 26      704  135. 14.6
# 27 27      367  137. 15.9
# 28 28      290  149. 10.9
# 29 29      366  136. 14.4
# 30 30      254  147. 11.2


