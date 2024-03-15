
library(tidyverse)
library(sweetsoursong)
library(spatstat.random)


# number of plants:
np <- 100

xy <- tibble(x = runif(np, 0, 10), y = runif(np, 0, 10))


ggplot(xy, aes(x, y)) +
    geom_point() +
    coord_equal(xlim = c(0, 10), ylim = c(0, 10))


z <- rStrauss(100, 1, 0.1)
z$n

tibble(x = z$x, y = z$y) |>
    ggplot(aes(x, y)) +
    geom_point() +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))

zz <- rMatClust(100, 0.01, 25)

tibble(x = zz$x, y = zz$y) |>
    ggplot(aes(x, y)) +
    geom_point() +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))



