#' Transition rates of the per-plant model
#'
#' Rates of the eight transitions in Table 1 for a plant with `y`
#' yeast-dominated, `b` bacteria-dominated, and `e` uncolonized flowers and
#' `p` pollinators. Arguments are recycled.
#'
#' @param y,b,p Numbers of yeast-dominated flowers, bacteria-dominated
#'   flowers, and pollinators.
#' @param e Number of uncolonized flowers. Defaults to `n - y - b`.
#' @param p0r,pyr,pbr Regional pools: mean numbers of pollinators per plant
#'   whose last visit was to an uncolonized, yeast-dominated, or
#'   bacteria-dominated flower, multiplied by `n`.
#' @param pars An [ss_params()] object.
#'
#' @return A matrix with one row per state and columns `p_minus`, `y_plus`,
#'   `b_plus`, `p_plus`, `py_plus`, `pb_plus`, `y_minus`, `b_minus`.
#' @export
flower_rates <- function(y, b, p, p0r, pyr, pbr, pars = ss_params(),
                         e = pars$n - y - b) {
    n <- pars$n; c <- pars$c; m <- pars$m; m_b <- pars$m_b
    e_y <- pars$e_y; e_b <- pars$e_b
    cbind(
        p_minus = c * p * (m * e / n + m * y / n + m_b * b / n),
        y_plus  = c * p * (1 - m) * e_y * y / n * e / n,
        b_plus  = c * p * (1 - m_b) * e_b * b / n * e / n +
            pars$c_b0 * b * e / n,
        p_plus  = c * (m * p0r / n + m * (1 - e_y * e / n) * pyr / n +
                           m_b * (1 - e_b * e / n) * pbr / n),
        py_plus = c * m * e_y * pyr / n * e / n,
        pb_plus = c * m_b * e_b * pbr / n * e / n,
        y_minus = pars$d * y,
        b_minus = pars$d * b
    )
}

#' Mean pollinators on a plant given its yeast-dominated flowers
#'
#' Balance of pollinator immigration and emigration in the reduced model:
#' `(m_b * pbr + m * pyr) / (m_b * (n - y) + m * y)`.
#'
#' @param y Number of yeast-dominated flowers (vectorized).
#' @param pyr,pbr Regional pools.
#' @param pars An [ss_params()] object.
#' @export
p_mean_given_y <- function(y, pyr, pbr, pars = ss_params()) {
    with(pars, (m_b * pbr + m * pyr) / (m_b * (n - y) + m * y))
}

#' Rates of the reduced model (no uncolonized flowers)
#'
#' A senescing flower is replaced immediately, by yeast or bacteria in
#' proportion to the rates at which an uncolonized flower would be colonized.
#' These are Chris Klausmeier's `fR` rates, with `B = n - Y`.
#'
#' @param y,p Numbers of yeast-dominated flowers and pollinators.
#' @param pyr,pbr Regional pools.
#' @param pars An [ss_params()] object.
#' @return A matrix with columns `p_minus`, `y_plus` (bacteria flower dies,
#'   local yeast fills), `b_plus`, `py_plus` (bacteria flower dies, immigrant
#'   yeast fills), `pb_plus`, `p_plus`.
#' @export
reduced_rates <- function(y, p, pyr, pbr, pars = ss_params()) {
    n <- pars$n; c <- pars$c; m <- pars$m; m_b <- pars$m_b
    b <- n - y
    # colonization rates of one uncolonized flower, divided by (e / n)
    ay_loc <- c * p * (1 - m) * pars$e_y * y / n
    ab_loc <- c * p * (1 - m_b) * pars$e_b * b / n + pars$c_b0 * b
    ay_imm <- c * m * pars$e_y * pyr / n
    ab_imm <- c * m_b * pars$e_b * pbr / n
    s <- ay_loc + ab_loc + ay_imm + ab_imm
    cbind(
        p_minus = c * p * (m * y / n + m_b * b / n),
        y_plus  = pars$d * b * ay_loc / s,
        b_plus  = pars$d * y * ab_loc / s,
        py_plus = pars$d * b * ay_imm / s,
        pb_plus = pars$d * y * ab_imm / s,
        p_plus  = c * (m * pyr / n + m_b * pbr / n)
    )
}
