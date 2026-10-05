#' Set small values to zero
#'
#' Equivalent to Wolfram Language `Chop`: values with absolute value below
#' `tol` become 0. Chris Klausmeier's Mathematica code applies `Chop` at
#' several steps; whether this package does the same is controlled by the
#' option `sweetsoursong.chop` (default `TRUE`, which reproduces the
#' Mathematica results).
#'
#' @param x Numeric vector or matrix.
#' @param tol Threshold.
#' @return `x` with small entries set to 0.
#' @export
chop <- function(x, tol = 1e-10) {
    x[abs(x) < tol] <- 0
    x
}

use_chop <- function() isTRUE(getOption("sweetsoursong.chop", TRUE))

# Chop only when the package option is on.
maybe_chop <- function(x) if (use_chop()) chop(x) else x
