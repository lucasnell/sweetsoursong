
#'
#' Preamble for use in all scripts.
#'

suppressPackageStartupMessages({
    library(sweetsoursong)
    library(tidyverse)
    library(patchwork)
})


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

#'
#' I add this to ggplot objects when I want to use them inside a figure
#' I'm stitching together in Adobe Illustrator, where I'll add all the titles
#' and annotations.
#'
subpanel_theme <- theme(plot.title = element_blank(),
                        legend.position = "none",
                        axis.title.y = element_blank(),
                        axis.title.x = element_blank(),
                        strip.text = element_blank())
