

library(sweetsoursong)
library(tidyverse)
library(patchwork)
library(RcppParallel)
setThreadOptions(numThreads = max(defaultNumThreads() - 2L, 1L))

if (file.exists(".Rprofile")) source(".Rprofile")

stoch_sims_file <- "_data/stoch-sims.rds"
open_stoch_sims_file <- "_data/stoch-open-sims.rds"


outcome_pal <- c("coexist" = "#008B00",
                 "yeast only" = "#FFCC33",
                 "bacteria only" = "#333399",
                 "extinct" = "gray60")
outcome_shapes <- c("coexist" = 19,
                    "yeast only" = 19,
                    "bacteria only" = 19,
                    "extinct" = 4)

spp_pal <- c(yeast = outcome_pal[["yeast only"]],
             bacteria = outcome_pal[["bacteria only"]],
             pollinators = "gray60")



#'
#' Simulations to compare dynamics across values of u and yeast dispersal
#' where one plant has lots of yeast and little bacteria, and another
#' plant is the opposite.
#'
# closed_sims <-
tibble(.u = c(0, 1, 1, 4), .d_yp = c(rep(1.4, 2), rep(1.9, 2)),
       pars = map(1:4, \(i) crossing()))


# |>
#     pmap_dfr(\(.u, .d_yp) {
#         high0 <- 0.5
#         low0 <- 0.02
#         plant_metacomm(np = 2, u = .u, d_yp = .d_yp,
#                        Y0 = c(high0, low0), B0 = c(low0, high0), max_t = 250) |>
#             mutate(u = .u, d_yp = .d_yp) |>
#             select(u, d_yp, everything())
#     }) |>
#     pivot_longer(Y:P, names_to = "type", values_to = "density") |>
#     mutate(type = factor(type, levels = c("Y", "B", "P"),
#                          labels = c("yeast", "bacteria", "pollinators")),
#            u = factor(u),
#            d_yp = factor(d_yp))






