
#' Create distance matrix from x and y coordinates
#'
#' @param x Numeric vector of locations in the x dimension.
#'     Must be at least 2 elements in length and must be the same length
#'     as argument `y`.
#'     Can be negative or positive, but cannot have infinite or missing values.
#'     Can also be a data frame with column names `"x"` and `"y"`.
#' @param y Numeric vector of locations in the y dimension.
#'     Must be at least 2 elements in length and must be the same length
#'     as argument `x`.
#'     Can be negative or positive, but cannot have infinite or missing values.
#'     Nothing should be passed to this argument if `x` is a data frame.
#'
#' @return A symmetrical numeric matrix with Euclidean distances between points.
#'
#' @export
#'
make_dist_mat <- function(x, y) {
    if (is.data.frame(x)) {
        stopifnot(is.data.frame(x) && all(c("x", "y") %in% colnames(x)))
        stopifnot(is.data.frame(x) && missing(y))
        y <- x$y
        x <- x$x
    }
    stopifnot(is.numeric(x) && is.numeric(y))
    stopifnot(length(x) == length(y))
    stopifnot(length(x) >= 2)
    stopifnot(all(!is.na(x)) && all(!is.na(y)))
    stopifnot(all(is.finite(x)) && all(is.finite(y)))
    dm <- make_dist_mat_rcpp(x, y)
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
    sw <- make_spat_wts_rcpp(dm, m)
    return(sw)
}




#

#' Create variance-covariance matrix from distance matrix
#'
#' @param dm Square, numeric distance matrix. Cannot have negative values.
#' @param q Single number indicating exponential distance decay constant.
#'     Must be >= 0.
#' @param sigma Single number or numeric vector of the same length as the
#'     number of rows in `dm`, indicating the standard deviations for
#'     all patches (if a single numeric) or for each patch individually
#'     (if a numeric vector).
#'     Cannot have negative values.
#'
#' @return A square, numeric variance-covariance matrix.
#'
#' @export
#'
make_vcv_mat <- function(dm, q, sigma) {
    stopifnot(is.matrix(dm) && is.numeric(dm) && nrow(dm) == ncol(dm))
    stopifnot(all(dm >= 0))
    stopifnot(is.numeric(q) && length(q) == 1 && q >= 0)
    stopifnot(is.numeric(sigma) && all(sigma >= 0))
    stopifnot(length(sigma) == 1 | length(sigma) == nrow(dm))
    if (length(sigma) == 1) sigma <- rep(sigma, nrow(dm))
    if (is.infinite(q)) {
        vcvm <- diag(nrow(dm))
    } else vcvm <- exp(-q * dm)
    diag(vcvm) <- sigma^2
    return(vcvm)
}

