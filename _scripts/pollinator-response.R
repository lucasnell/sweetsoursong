#'
#' How do pollinators respond to flower patch size?
#'
#' This dataset spans 9 years but cannot be used for phenology because they
#' didn't sample late into the flowering season.
#'
#' These analyses didn't work out as interesting because the data are just too
#' noisy. It probably has to do with how they were collected since touching
#' the flowers can cause their stigmas to close.
#'
#' I also tried analyzing yeast CFU as a proxy for pollinator activity, but
#' that's noisy too.
#'

library(vegan)
library(lme4)
library(tidyverse)

theme_set(theme_classic())

logit <- function(p) log(p/(1-p))
inv_logit <- function(x) 1 / (1 + exp(-x))





date_df <- paste0("~/Stanford_Drive/_sweetandsour/BIO47-data-2024-03-11/",
                  "2023_main_tables/week_dates.csv") |>
    read_csv(col_types = "icccc") |>
    rename(year = Year) |>
    mutate(across(start.date:end.date,
                  \(x) {
                      str_split(x, "\\s+") |>
                          map_chr(\(x) x[[1]]) |>
                          as.Date(format = "%m/%d/%y")
                  })) |>
    filter(grepl("Week", period.name)) |>
    mutate(week = pmap(list(period.name, start.date, end.date),
                       \(x, s_d, e_d) {
                           if (grepl("and", x)) {
                               w <- str_remove(x, "Weeks") |>
                                   str_split("and") |>
                                   getElement(1) |>
                                   as.integer()
                               m_d <- median(c(s_d, e_d))
                               tibble(week = w,
                                      start = c(s_d, m_d),
                                      end = c(m_d, e_d))
                           } else {
                               w <- as.integer(gsub("Week", "", x))
                               tibble(week = w, start = s_d, end = e_d)
                           }
                       })) |>
    unnest(week) |>
    select(year, week, start, end)


#' Get middle of week and account for missing weeks in `date_df`.
#' This is already vectorized and should not be used inside a `map` or
#' `apply` function.
mow_getter <- function(year, week) {
    weeks <- date_df$week
    years <- date_df$year
    starts <- date_df$start
    ends <- date_df$end
    # Middles of the weeks:
    mows <- map2_vec(year, week, \(yr, wk) {
        idx <- which(years == yr & weeks == wk)
        if (length(idx) == 0) return(NA_Date_)
        if (length(idx) > 1) return(NA_Date_)
        mow <- mean(c(starts[idx], ends[idx]))
        return(mow)
    })
    # fill in those that don't have direct matches:
    mows[is.na(mows)] <- map_vec(which(is.na(mows)),
                                 \(i) {
                                     yr <- year[i]
                                     wk <- week[i]
                                     yw_diffs <- wk - weeks[years == yr]
                                     mdi <- which(abs(yw_diffs) == min(abs(yw_diffs)))
                                     s_d <- starts[years == yr][mdi] + yw_diffs[mdi] * 7
                                     e_d <- ends[years == yr][mdi] + yw_diffs[mdi] * 7
                                     return(mean(c(s_d, e_d)))
                                 })
    return(mows)
}


flower_df <- paste0("~/Stanford_Drive/_sweetandsour/BIO47-data-2024-03-11/",
                    "2023_main_tables/flower_counts.csv") |>
    read_csv(col_types = "iiidddd") |>
    rename_with(\(x) tolower(x) |>
                    str_remove_all("_stigmas") |>
                    str_replace_all("num", "n") |>
                    str_replace_all("fraction", "p")) |>
    mutate(date = mow_getter(year, week),
           doy = yday(date),
           p_closed = ifelse(n_flowers == 0, NA, p_closed)) |>
    arrange(year, plant, week) |>
    group_by(year, plant) |>
    mutate(n_flowers_lag = lag(n_flowers),
           n_flowers_lag2 = lag(n_flowers, 2),
           n_flowers_lag3 = lag(n_flowers, 3)) |>
    ungroup() |>
    # Plotting this year shows that something clearly changed with the methods:
    filter(year != 2023) |>
    mutate(across(year:plant, factor))


cfu_df <- paste0("~/Stanford_Drive/_sweetandsour/BIO47-data-2024-03-11/",
                 "2023_main_tables/flowers.csv") |>
    read_csv(col_types = cols()) |>
    rename_with(tolower) |>
    rename(yeast_cfu_ul = fungal_cfu_per_ul) |>
    filter(!bagged, !caged, !is.na(yeast_cfu_ul), !is.na(week)) |>
    select(plant:flower, stigma, nectar_ph, nectar_ul, yeast_cfu_ul,
           age_category, age_days) |>
    mutate(n_flowers_lag =
               pmap_dbl(list(plant, year, week),
                        \(pl, yr, wk) {
                            lgl <- flower_df$plant == pl &
                                flower_df$year == yr &
                                flower_df$week == wk
                            out <- flower_df$n_flowers_lag[lgl]
                            if (length(out) == 0) return(NA)
                            return(out)
                        })) |>
    filter(!is.na(n_flowers_lag)) |>
    mutate(id = interaction(year, plant, drop = TRUE, sep = "_"),
           log_n_flowers_lag = log1p(n_flowers_lag),
           has_yeast = yeast_cfu_ul > 0)


cfu_df |>
    ggplot(aes(log_n_flowers_lag, log1p(yeast_cfu_ul))) +
    geom_point(shape = 1) +
    facet_wrap(~year) +
    xlab("log(Number of flowers in previous week + 1)") +
    ylab("log(Yeast CFU / µL + 1)") +
    theme(strip.background = element_rect(color = NA),
          strip.text = element_text(size = 11, face = "bold"))

cfu_df |>
    ggplot(aes(log_n_flowers_lag, log1p(yeast_cfu_ul))) +
    geom_point(shape = 1) +
    xlab("log(Number of flowers in previous week + 1)") +
    ylab("log(Yeast CFU / µL + 1)")
cfu_df |>
    ggplot(aes(log_n_flowers_lag, has_yeast)) +
    geom_point(shape = 1) +
    stat_summary(fun.data = "mean_cl_boot", color = "dodgerblue",
                 linewidth = 2, size = 1, shape = 3) +
    xlab("log(Number of flowers in previous week + 1)") +
    ylab("Yeast presence")



cfu_mod <- lmer(log1p(yeast_cfu_ul) ~ log_n_flowers_lag + I(log_n_flowers_lag^2) +
                    (1 | id), cfu_df,
                REML = FALSE)
cfu_mod |> AIC()
# best AIC = 1095.54 for Y ~ log_n_flowers_lag + I(log_n_flowers_lag^2) + (1 | id)
cfu_mod |> summary()
cfu_mod |> car::Anova()

pres_mod <- glmer(has_yeast ~ log_n_flowers_lag + I(log_n_flowers_lag^2) +
                    (1 | id),
                  family = binomial, cfu_df)
pres_mod |> AIC()
# best AIC = 283.2563 for Y ~ log_n_flowers_lag + I(log_n_flowers_lag^2) + (1 | id)
pres_mod |> summary()
pres_mod |> car::Anova()















flower_df |>
    filter(!is.na(n_flowers_lag), !is.na(p_closed)) |>
    ggplot(aes((n_flowers_lag), (p_closed * n_flowers))) +
    geom_point() +
    scale_y_continuous("Number of stigma closed", trans = "log1p") +
    scale_x_continuous("Number of flowers in previous week", trans = "log1p")
flower_df |>
    filter(!is.na(n_flowers), !is.na(p_closed)) |>
    ggplot(aes(n_flowers, (p_closed))) +
    geom_point() +
    scale_y_continuous("Proportion of stigma closed", trans = "logit") +
    scale_x_continuous("Number of flowers in current week", trans = "log1p")
flower_df |>
    filter(!is.na(n_flowers_lag), !is.na(p_closed)) |>
    ggplot(aes(n_flowers_lag, (p_closed))) +
    geom_point() +
    scale_y_continuous("Proportion of stigma closed", trans = "logit") +
    scale_x_continuous("Number of flowers in previous week", trans = "log1p")





# Flower numbers rounded to nearest integer:
int_flower_df <- flower_df |>
    filter(!is.na(n_flowers_lag), !is.na(p_closed)) |>
    mutate(across(n_closed:n_open, \(x) round(x)),
           log_n_flowers = log1p(n_flowers),
           log_n_flowers_lag = log1p(n_flowers_lag),
           log_n_flowers_lag2 = log1p(n_flowers_lag2),
           log_n_flowers_lag3 = log1p(n_flowers_lag3),
           id = interaction(year, plant, drop = TRUE, sep = "_"))

mod <- glmer(cbind(n_closed, n_open) ~ log_n_flowers_lag +
                 (log_n_flowers_lag | id),
             family = binomial, int_flower_df)

mod |> AIC()
# best AIC = 3005.317 for Y ~ log_n_flowers_lag + (log_n_flowers_lag | id)
mod |> summary()

mod0 <- glmer(cbind(n_closed, n_open) ~ log_n_flowers +
                  (1 | id),
             family = binomial, int_flower_df)

mod0 |> AIC()
# best AIC = 2992.768 for Y ~ log_n_flowers + (1 | id)
mod0 |> summary()


mod2 <- glmer(cbind(n_closed, n_open) ~ log_n_flowers_lag2 +
                 (1 | id),
             family = binomial,
             filter(int_flower_df, !is.na(log_n_flowers_lag2)))

mod2 |> AIC()
# best AIC = 2544.193 for Y ~ log_n_flowers_lag2 + (1 | id)


mod3 <- glmer(cbind(n_closed, n_open) ~ log_n_flowers_lag3 +
                 (1 | id),
             family = binomial,
             filter(int_flower_df, !is.na(log_n_flowers_lag3)))






pred_df <- tibble(log_n_flowers_lag = int_flower_df$log_n_flowers_lag |>
           (\(x) seq(min(x), max(x), length.out = 101))()) |>
    (\(nd) {
        z <- suppressWarnings(predict(mod, re.form = NA, newdata = nd,
                                      type = "response", se.fit = TRUE))
        mutate(nd, p_closed = z$fit, p_closed_se = z$se.fit)
    })() |>
    mutate(n_flowers_lag = exp(log_n_flowers_lag) - 1)

pred_df2 <- tibble(log_n_flowers_lag2 = int_flower_df$log_n_flowers_lag2 |>
           (\(x) seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 101))()) |>
    (\(nd) {
        z <- suppressWarnings(predict(mod2, re.form = NA, newdata = nd,
                                      type = "response", se.fit = TRUE))
        mutate(nd, p_closed = z$fit, p_closed_se = z$se.fit)
    })() |>
    mutate(n_flowers_lag2 = exp(log_n_flowers_lag2) - 1)




int_flower_df |>
    ggplot(aes(log_n_flowers_lag, p_closed)) +
    geom_ribbon(aes(x = log_n_flowers_lag, ymin = p_closed - p_closed_se,
                    ymax = p_closed + p_closed_se),
                data = pred_df, fill = "gray80") +
    geom_line(data = pred_df, color = "dodgerblue", linewidth = 1) +
    geom_point(aes(color = year)) +
    # facet_wrap(~ year, nrow = 3) +
    scale_color_viridis_d(guide = "none") +
    ylab("Proportion stigma closed") +
    xlab("log(Number of flowers in previous week + 1)")


int_flower_df |>
    filter(!is.na(n_flowers_lag2)) |>
    ggplot(aes(log_n_flowers_lag2, p_closed)) +
    geom_ribbon(aes(x = log_n_flowers_lag2, ymin = p_closed - p_closed_se,
                    ymax = p_closed + p_closed_se),
                data = pred_df2, fill = "gray80") +
    geom_line(data = pred_df2, color = "dodgerblue", linewidth = 1) +
    geom_point(aes(color = year)) +
    # facet_wrap(~ year, nrow = 3) +
    scale_color_viridis_d(guide = "none") +
    ylab("Proportion stigma closed") +
    xlab("log(Number of flowers 2 weeks ago + 1)")






mod0 <- glm(cbind(n_closed, n_open) ~ log_n_flowers_lag, family = binomial, int_flower_df)
mod1 <- glm(cbind(n_closed, n_open) ~ log_n_flowers_lag + I(log_n_flowers_lag^2), family = binomial, int_flower_df)

# Likelihood ratio test:
DF <- attr(logLik(mod1), "df") - attr(logLik(mod0), "df")
LR <- -2 * (as.numeric(logLik(mod0)) - as.numeric(logLik(mod1)))
pchisq(LR, DF, lower.tail = FALSE)


mos_mod <- MOStest(int_flower_df$log_n_flowers_lag,
                   cbind(int_flower_df$n_closed, int_flower_df$n_open),
                   family = binomial(link = "logit"))
# plot(mos_mod)
mos_mod



