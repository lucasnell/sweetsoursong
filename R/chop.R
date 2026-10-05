#' Set small values to zero
#'
#' Equivalent to Wolfram Language `Chop`: values with absolute value below
#' `tol` become 0. Chris Klausmeier's Mathematica code applies `Chop` at
#' several steps. Set the option `sweetsoursong.chop = TRUE` to do the same
#' and reproduce his results exactly. The default, `FALSE`, does not, because
#' Chop zeroes small probabilities that matter for invasion criteria: in the
#' reduced model with `eps <= 1e-7` it sets the yeast invasion criterion to
#' 0, and in the full model it shifts the invasion thresholds by up to 0.005
#' in `pr`. Without Chop, negative round-off is set to 0.
#'
#' @param x Numeric vector or matrix.
#' @param tol Threshold.
#' @return `x` with small entries set to 0.
#' @export
chop <- function(x, tol = 1e-10) {
    x[abs(x) < tol] <- 0
    x
}

use_chop <- function() isTRUE(getOption("sweetsoursong.chop", FALSE))

# Chop only when the package option is on.
maybe_chop <- function(x) if (use_chop()) chop(x) else x
