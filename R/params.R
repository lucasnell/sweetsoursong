#' Model parameters
#'
#' Parameter values for the per-plant model. Defaults are `SetParameters` in
#' Chris Klausmeier's `exploration_unified_10x.nb`, which produce the
#' manuscript figures.
#'
#' @param n Number of flowers per plant (`N` in the manuscript).
#' @param pmax Maximum number of pollinators per plant in numerical solutions.
#' @param c Rate at which each pollinator visits flowers.
#' @param d Flower senescence rate.
#' @param m Probability that a pollinator leaves the plant after visiting an
#'   uncolonized or yeast-dominated flower.
#' @param m_b Probability that a pollinator leaves the plant after visiting a
#'   bacteria-dominated flower.
#' @param e_y,e_b Probabilities that a pollinator arriving from a yeast- or
#'   bacteria-dominated flower establishes that microbe in an uncolonized
#'   flower.
#' @param c_b0 Rate at which bacteria colonize uncolonized flowers without
#'   pollinators, per bacteria-dominated flower.
#' @param eps Size of the invader's regional pool for invasion criteria.
#'
#' @return A named list of class `ss_params`.
#' @export
ss_params <- function(n = 50L, pmax = 12L, c = 500, d = 0.1, m = 0.01,
                      m_b = 0.05, e_y = 1, e_b = 0.5, c_b0 = 5, eps = 1e-5) {
    pars <- list(n = as.integer(n), pmax = as.integer(pmax), c = c, d = d,
                 m = m, m_b = m_b, e_y = e_y, e_b = e_b, c_b0 = c_b0,
                 eps = eps)
    check_params(pars)
    class(pars) <- "ss_params"
    pars
}

#' Change parameter values
#'
#' @param pars An `ss_params` object.
#' @param ... Named values to replace.
#' @return An `ss_params` object.
#' @export
update_params <- function(pars, ...) {
    new <- list(...)
    bad <- setdiff(names(new), names(pars))
    if (length(bad) > 0) stop("unknown parameters: ", paste(bad, collapse = ", "))
    pars[names(new)] <- new
    do.call(ss_params, unclass(pars))
}

check_params <- function(pars) {
    with(pars, {
        stopifnot(
            length(n) == 1, n >= 1,
            length(pmax) == 1, pmax >= 1,
            c > 0, d > 0, c_b0 >= 0, eps > 0,
            m >= 0, m <= 1, m_b >= 0, m_b <= 1,
            e_y >= 0, e_y <= 1, e_b >= 0, e_b <= 1
        )
    })
    invisible(TRUE)
}

#' @export
print.ss_params <- function(x, ...) {
    cat("sweetsoursong parameters\n")
    v <- unlist(unclass(x))
    print(v, ...)
    invisible(x)
}

#' Critical pollinator number
#'
#' With pollinator numbers fixed and no microbe immigration, yeast excludes
#' bacteria on a plant when `P > p_crit` and bacteria exclude yeast when
#' `P < p_crit`.
#'
#' @inheritParams e_yp
#' @return `c_b0 * n / (c * (e_y - e_b))`.
#' @export
p_crit <- function(pars = ss_params()) {
    with(pars, c_b0 * n / (c * (e_y - e_b)))
}
