#' Stationary distribution of a birth-death chain
#'
#' Product formula `pi(y) = pi(y - 1) * up(y - 1) / down(y)`, rescaled when
#' values exceed 1e100 (Chris Klausmeier's `CTMCClosureStationaryWeights`).
#'
#' @param up,down Birth and death rates for states 0, ..., n.
#' @return Probabilities for states 0, ..., n.
#' @export
birth_death_stationary <- function(up, down) {
    k <- length(up)
    stopifnot(length(down) == k)
    w <- numeric(k)
    w[1] <- 1
    for (i in seq_len(k - 1) + 1L) {
        w[i] <- if (down[i] > 0) w[i - 1] * up[i - 1] / down[i] else 0
        s <- max(w)
        if (s > 1e100) w <- w / s
    }
    if (sum(w) > 0) w / sum(w) else w
}

#' Stationary distribution of yeast-dominated flowers on one plant
#'
#' Reduced model with the vacancy closure, for fixed regional pools
#' (Chris Klausmeier's `GetYDistribution` with `SetCTMCClosureTM`).
#'
#' @param pyr,pbr Regional pools.
#' @param pars An [ss_params()] object.
#' @return A numeric vector of probabilities for y = 0, ..., n, with the
#'   hazards from [vacancy_hazards()] as attribute `hazards`.
#' @export
y_distribution <- function(pyr, pbr, pars = ss_params()) {
    hz <- vacancy_hazards(pyr, pbr, pars)
    pi_y <- birth_death_stationary(hz$up, hz$down)
    attr(pi_y, "hazards") <- hz
    pi_y
}

#' Joint stationary distribution of yeast-dominated flowers and pollinators
#'
#' `pi(y) * Poisson(p; E[P | y])` for p = 0, ..., `pmax`, normalized over the
#' truncated grid (Chris Klausmeier's `GetPDF`).
#'
#' @inheritParams y_distribution
#' @param pi_y Optional output of [y_distribution()] for the same pools.
#' @return An (n + 1) x (pmax + 1) matrix; rows are y = 0, ..., n and columns
#'   p = 0, ..., pmax.
#' @export
joint_distribution <- function(pyr, pbr, pars = ss_params(),
                               pi_y = y_distribution(pyr, pbr, pars)) {
    n <- pars$n
    lambda <- p_mean_given_y(0:n, pyr, pbr, pars)
    pois <- t(vapply(lambda, function(l) stats::dpois(0:pars$pmax, l),
                     numeric(pars$pmax + 1)))
    joint <- as.numeric(pi_y) * pois
    joint <- joint / sum(joint)
    dimnames(joint) <- list(y = 0:n, p = 0:pars$pmax)
    joint
}

#' Means of a joint distribution over (y, p)
#'
#' @param joint Output of [joint_distribution()].
#' @return Named vector with means of y and p.
#' @export
joint_means <- function(joint) {
    y <- as.numeric(rownames(joint))
    p <- as.numeric(colnames(joint))
    tot <- sum(joint)
    c(y = sum(rowSums(joint) * y) / tot, p = sum(colSums(joint) * p) / tot)
}
