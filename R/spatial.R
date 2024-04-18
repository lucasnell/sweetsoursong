
#' Create distance matrix from x and y coordinates
#'
#' @param x Numeric vector of locations in the x dimension.
#'     Must be at least 2 elements in length and must be the same length
#'     as argument `y`.
#'     Can be negative or positive, but cannot have infinite or missing values.
#' @param y Numeric vector of locations in the y dimension.
#'     Must be at least 2 elements in length and must be the same length
#'     as argument `x`.
#'     Can be negative or positive, but cannot have infinite or missing values.
#'
#' @return A symmetrical numeric matrix with Euclidean distances between points.
#'
#' @export
#'
make_dist_mat <- function(x, y) {
    stopifnot(is.numeric(x) && is.numeric(y))
    stopifnot(length(x) == length(y))
    stopifnot(length(x) >= 2)
    stopifnot(all(!is.na(x)) && all(!is.na(y)))
    stopifnot(all(is.finite(x)) && all(is.finite(y)))
    n <- length(x)
    dm <- matrix(0, n, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            dm[i,j] <- dm[j,i] <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
        }
    }
    return(dm)
}


#' Create a matrix of spatial weights from distance matrix
#'
#' @param dm Numeric distance matrix. Must be symmetrical with zeros on
#'     the diagonal and only contain non-negative values.
#'     Missing values are not allowed.
#' @param m A single numeric indicating the exponent for spatial weighting.
#'     Spatial weighting is `d^m` for each item `d` in `dm`.
#'     The default is `2`.
#'
#' @return A symmetrical numeric matrix with spatial weighting between points.
#'
#' @export
#'
make_spat_wts <- function(dm, m = 2) {
    stopifnot(is.numeric(m) && length(m) == 1 && is.finite(m))
    stopifnot(is.matrix(dm) && is.numeric(dm) && isSymmetric(dm) && nrow(dm) > 1)
    stopifnot(all(!is.na(dm)) && all(is.finite(dm)) && all(diag(dm) == 0))
    n <- nrow(dm)
    sw <- matrix(0, n, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            sw[i,j] <- sw[j,i] <- 1 / (dm[i,j]^m)
        }
    }
    return(sw)
}

