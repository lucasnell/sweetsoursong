
library(tidyverse)
library(sweetsoursong)
library(spatstat.random)





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

zz <- rMatClust(kappa = 20, scale = 0.05, mu = 5)
zz$n

tibble(x = zz$x, y = zz$y) |>
    ggplot(aes(x, y)) +
    geom_point() +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))





# create distance matrix:
make_dist_mat <- function(x, y) {
    stopifnot(length(x) == length(y))
    stopifnot(length(x) >= 2)
    n <- length(x)
    dm <- matrix(0, n, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            dm[i,j] <- dm[j,i] <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
        }
    }
    return(dm)
}

# from distance to spatial weights:
make_spat_wts <- function(dm, m = 1) {
    stopifnot(is.matrix(dm) && is.numeric(dm) && isSymmetric(dm) && nrow(dm) > 1)
    n <- nrow(dm)
    sw <- matrix(0, n, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            sw[i,j] <- sw[j,i] <- 1 / (dm[i,j]^m)
        }
    }
    return(sw)
}
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
X <- MASS::mvrnorm(n = 1, mu = rep(1, nrow(vcv_mat)), var_cov_mat(dist_mat, 1, 0.5))
# Calculate Moran's I (spatial autocorrelation measure)
autocor(X, ww, "moran")

# q = 10 should reduce Moran's I
X <- MASS::mvrnorm(n = 1, mu = rep(1, nrow(vcv_mat)), var_cov_mat(dist_mat, 10, 0.5))
autocor(X, ww, "moran")


