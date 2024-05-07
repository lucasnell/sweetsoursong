

suppressPackageStartupMessages({
    library(MASS)  # mvrnorm
    library(terra) # autocor
    library(spatstat.random) # rMatClust
    library(sweetsoursong)
    library(tidyverse)
    library(patchwork)
    library(gganimate)
    library(parallel)
})


make_Phi <- function(dm, w) {
    stopifnot(is.matrix(dm) && is.numeric(dm))
    stopifnot(nrow(dm) == ncol(dm))
    stopifnot(is.numeric(w) && length(w) == 1)
    n_plants <- nrow(dm)
    Phi <- matrix(0, n_plants, n_plants)
    for (j in 1:n_plants) {
        col_sum = 0
        for (i in 1:n_plants) {
            if (i == j) {
                Phi[i,j] = 1
            } else {
                Phi[i,j] = exp(-w * dm[i, j])
            }
            col_sum = col_sum + Phi[i,j]
        }
        for (i in 1:n_plants) Phi[i,j] = Phi[i,j] / col_sum
    }
    return(Phi)
}
make_Phi_disk <- function(dm, a) {
    stopifnot(is.matrix(dm) && is.numeric(dm))
    stopifnot(nrow(dm) == ncol(dm))
    stopifnot(is.numeric(a) && length(a) == 1)
    n_plants <- nrow(dm)
    Phi <- matrix(0, n_plants, n_plants)
    for (j in 1:n_plants) {
        col_sum = 0
        for (i in 1:n_plants) {
            if (i == j) {
                Phi[i,j] = 1
            } else if (dm[i, j] > a) {
                Phi[i,j] = 0
            } else {
                da = dm[i, j] / a
                Phi[i,j] = (2 / pi) * { acos(da) - sqrt(da * { 1 - da^2 }) }
            }
            col_sum = col_sum + Phi[i,j]
        }
        for (i in 1:n_plants) Phi[i,j] = Phi[i,j] / col_sum
    }
    return(Phi)
}
getP <- function(w, B) {
    wts <- numeric(nrow(dm))
    for(i in 1:nrow(dm)) {
        for (j in 1:ncol(dm)) {
            wts[i] <- wts[i] + B[j] * exp(-w * dm[i,j])
        }
    }
    wts <- wts / sum(wts)
    return(wts)
}
getP_disk <- function(a, B) {
    wts <- numeric(nrow(dm))
    for(i in 1:nrow(dm)) {
        d <- dm[i,]
        d[d > a] <- 0
        d[d <= a] <- d[d <= a] / a
        d[d <= a] <- (2 / pi) * { acos(d[d <= a]) -
                sqrt(d[d <= a] * { 1 - (d[d <= a])^2 }) }
        for (j in 1:ncol(dm)) {
            wts[i] <- wts[i] + B[j] * d[j]
        }
    }
    wts <- wts / sum(wts)
    return(wts)
}


dm <- rMatClust(kappa = 25, scale = 0.01, mu = 5, nsim = 1) |>
    as.data.frame() |>
    slice_sample(n = 100) |>
    make_dist_mat()
stopifnot(nrow(dm) == 100)

a__ <- 0.01
w__ <- 1


Phi <- make_Phi(dm, w = w__)
Phi2 <- make_Phi_disk(dm, a = a__)
# B <- rep(1, 100)
B <- runif(100)
# B <- mvrnorm(1, rep(100, 100), make_vcv_mat(dm, 1, q = 10))

P <- getP(w = w__, B)
P_disk <- getP_disk(a = a__, B)
P2 <- Phi %*% cbind(B)
P3 <- Phi2 %*% cbind(B)
# par(mfrow = c(1, 4)); hist(P, breaks = 30); hist(P_disk, breaks = 30); hist(P2, breaks = 30); hist(P3, breaks = 30)


for (i in 1:100) {
    P <- getP(w = 1, P)
    P_disk <- getP_disk(a = a__, P_disk)
    P2 <- Phi %*% P2
    P3 <- Phi2 %*% P3
}
# par(mfrow = c(1, 4)); hist(P, breaks = 30); hist(P_disk, breaks = 30); hist(P2, breaks = 30); hist(P3, breaks = 30)

dist_Pdiff <- tibble(pw_dist = map(2:nrow(dm),\(i) map_dbl(1:(i-1),\(j) dm[i,j])) |> do.call(what = c))
dist_Pdiff$P <-  map(2:nrow(dm),\(i) map_dbl(1:(i-1),\(j) abs(P[i] - P[j]))) |> do.call(what = c)
dist_Pdiff$P_disk <-  map(2:nrow(dm),\(i) map_dbl(1:(i-1),\(j) abs(P_disk[i] - P_disk[j]))) |> do.call(what = c)
dist_Pdiff$P2 <-  map(2:nrow(dm),\(i) map_dbl(1:(i-1),\(j) abs(P2[i] - P2[j]))) |> do.call(what = c)
dist_Pdiff$P3 <-  map(2:nrow(dm),\(i) map_dbl(1:(i-1),\(j) abs(P3[i] - P3[j]))) |> do.call(what = c)

dist_Pdiff |>
    mutate(across(P:P3, \(x) x / max(x))) |>
    pivot_longer(P:P3) |>
    ggplot(aes(pw_dist, value)) +
    geom_point(shape = 1, alpha = 0.5) +
    stat_smooth(method = "lm", formula = y ~ x) +
    facet_wrap(~ name, nrow = 2)



X <- cbind(runif(3))
X

M <- matrix(0, 3, 3)
M[lower.tri(M)] <- c(0.01, 0.1, 0.6)
M <- M + t(M) + diag(3)
V <- M / matrix(colSums(M), 3, 3, byrow=TRUE)

X; V %*% X; V

Z <- V %*% X
for (i in 1:1e6L) Z <- V %*% Z
Z
V
V |> rowSums()
eigen(V)
eigen(t(V))


mv_df <- tibble(t = 0:10e3L, X1 = 100, X2 = -100, X3 = 0)
for (i in 2:nrow(mv_df)) {
    ZZ <- c(mv_df$X1[[i-1]], mv_df$X2[[i-1]], mv_df$X3[[i-1]])
    XX <- ZZ + mvrnorm(1, rep(0, 3), M - diag(3) * 0.3)
    mv_df$X1[[i]] <- XX[1]
    mv_df$X2[[i]] <- XX[2]
    mv_df$X3[[i]] <- XX[3]
}; rm(i, XX, ZZ)
mv_df |>
    pivot_longer(X1:X3) |>
    ggplot(aes(t, value, color = name)) +
    geom_line()











make_Phi <- function(dm, w) {
    stopifnot(is.matrix(dm) && is.numeric(dm))
    stopifnot(nrow(dm) == ncol(dm))
    stopifnot(is.numeric(w) && length(w) == 1)
    n_plants <- nrow(dm)
    Phi <- matrix(0, n_plants, n_plants)
    for (j in 1:n_plants) {
        col_sum = 0
        for (i in 1:n_plants) {
            if (i == j) {
                Phi[i,j] = 1
            } else {
                Phi[i,j] = exp(-w * dm[i, j])
            }
            col_sum = col_sum + Phi[i,j]
        }
        for (i in 1:n_plants) Phi[i,j] = Phi[i,j] / col_sum
    }
    return(Phi)
}

make_Phi2 <- function(dm, a) {
    stopifnot(is.matrix(dm) && is.numeric(dm))
    stopifnot(nrow(dm) == ncol(dm))
    stopifnot(is.numeric(a) && length(a) == 1)
    f <- function(d, a) {
        if (d > a) return(0)
        (2 / pi) * (acos(d / a) - (d/a * (1 - (d/a)^2))^0.5)
    }
    n_plants = nrow(dm)
    Phi = matrix(0, n_plants, n_plants)
    for (j in 1:n_plants) {
        col_mean = 0
        for (i in 1:n_plants) {
            if (i == j) {
                Phi[i,j] = 1
            } else {
                Phi[i,j] = f(dm[i,j], a)
            }
            col_sum = col_sum + Phi[i,j]
        }
        for (i in 1:n_plants) Phi[i,j] = Phi[i,j] / col_sum
    }
    return(Phi)
}

np <- 3L
set.seed(92345678)
Y <- runif(np, 0, 150)
dm <- make_dist_mat(runif(np), runif(np))

Phi <- make_Phi(dm, w = 5)
Phi2 <- make_Phi2(dm, a = 0.25)

Phi %*% cbind(Y); cbind(Y)
sum(Phi %*% cbind(Y)); sum(Y)

Phi2 %*% cbind(Y); cbind(Y)


#


