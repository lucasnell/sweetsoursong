
library(MASS)  # mvrnorm
library(terra) # autocor
library(spatstat.random) # rMatClust
library(tidyverse)
library(sweetsoursong)





# number of plants:
np <- 100

xy <- tibble(x = runif(np, 0, 10), y = runif(np, 0, 10))


ggplot(xy, aes(x, y)) +
    geom_point() +
    coord_equal(xlim = c(0, 10), ylim = c(0, 10))


# z <- rStrauss(100, gamma = 0.9, R = 0.2)
# z$n
#
# tibble(x = z$x, y = z$y) |>
#     ggplot(aes(x, y)) +
#     geom_point() +
#     coord_equal(xlim = c(0, 1), ylim = c(0, 1))

zz <- rMatClust(kappa = 20, scale = 0.05, mu = 5, nsim = 1)
zz$n

tibble(x = zz$x, y = zz$y) |>
    ggplot(aes(x, y)) +
    geom_point() +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))




# to variance-covariance matrix
var_cov_mat <- function(dm, q, sigma) {
    stopifnot(is.numeric(q) && length(q) == 1)
    stopifnot(is.numeric(sigma) && length(sigma) == 1)
    vcvm <- sigma * exp(-q * dm)
    return(vcvm)
}

dist_mat <- make_dist_mat(zz$x, zz$y)
ww <- make_spat_wts(dist_mat)

# Generate spatially autocorrelated random variable with q = 1:
X <- mvrnorm(n = 1, mu = rep(1, nrow(dist_mat)), var_cov_mat(dist_mat, 1, 0.5))
# Calculate Moran's I (spatial autocorrelation measure)
autocor(X, ww, "moran")

# q = 10 should reduce Moran's I
X <- mvrnorm(n = 1, mu = rep(1, nrow(dist_mat)), var_cov_mat(dist_mat, 10, 0.5))
autocor(X, ww, "moran")


