#' Deterministic one-plant model
#'
#' Right-hand side of the reduced model (equation 2) for fixed regional
#' pools. With `immigration = FALSE`, pollinators move between plants but do
#' not carry microbes (Chris Klausmeier's `SetNoImmigrationEcoEvoModel`);
#' with `TRUE`, they do (`SetReducedEcoEvoModel`).
#'
#' @param y,p Yeast-dominated flowers and pollinators.
#' @param pyr,pbr Regional pools.
#' @param pars An [ss_params()] object.
#' @param immigration Whether arriving pollinators carry microbes.
#' @return A matrix with columns `dy` and `dp`.
#' @export
determ_rhs <- function(y, p, pyr, pbr, pars = ss_params(), immigration = TRUE) {
    r <- reduced_rates(y, p, pyr, pbr, pars)
    dy <- r[, "y_plus"] - r[, "b_plus"]
    if (immigration) dy <- dy + r[, "py_plus"] - r[, "pb_plus"]
    cbind(dy = dy, dp = r[, "p_plus"] - r[, "p_minus"])
}

#' Simulate the deterministic one-plant model
#'
#' @inheritParams determ_rhs
#' @param y0,p0 Initial values.
#' @param times Output times.
#' @return A data frame from [deSolve::ode()].
#' @export
determ_sim <- function(y0, p0, times, pyr, pbr, pars = ss_params(),
                       immigration = TRUE) {
    fn <- function(t, x, parms) {
        list(as.numeric(determ_rhs(x[1], x[2], pyr, pbr, pars, immigration)))
    }
    out <- deSolve::ode(c(y = y0, p = p0), times, fn, parms = NULL)
    as.data.frame(out)
}

#' Equilibria of the deterministic one-plant model
#'
#' Pollinators equilibrate at `p = p_mean_given_y(y)`; equilibria are the
#' roots in `[0, n]` of `dy` along that curve, found by a grid scan and
#' [stats::uniroot()]. Stability is from the eigenvalues of a numerical
#' Jacobian.
#'
#' @inheritParams determ_rhs
#' @param n_grid Number of grid points for the scan.
#' @return A data frame with `y`, `p`, `stable`, and the two eigenvalues
#'   (real parts `ev1_re`, `ev2_re`; imaginary parts `ev1_im`, `ev2_im`).
#' @export
determ_equilibria <- function(pyr, pbr, pars = ss_params(), immigration = TRUE,
                              n_grid = 20001L) {
    n <- pars$n
    g <- function(y) determ_rhs(y, p_mean_given_y(y, pyr, pbr, pars), pyr,
                                pbr, pars, immigration)[, "dy"]
    ys <- seq(0, n, length.out = n_grid)
    gs <- g(ys)
    roots <- ys[gs == 0]
    sc <- which(sign(gs[-1]) * sign(gs[-n_grid]) < 0)
    for (i in sc) {
        roots <- c(roots, stats::uniroot(g, ys[c(i, i + 1)], tol = 1e-13)$root)
    }
    roots <- sort(unique(roots))
    jac <- function(y, p) {
        h <- 1e-6
        f <- function(z) as.numeric(determ_rhs(z[1], z[2], pyr, pbr, pars, immigration))
        x <- c(y, p)
        J <- matrix(0, 2, 2)
        for (k in 1:2) {
            hk <- h * max(1, abs(x[k]))
            e <- c(0, 0); e[k] <- hk
            lo <- x - e
            # one-sided difference at the boundary y = 0
            if (lo[1] < 0) J[, k] <- (f(x + e) - f(x)) / hk
            else J[, k] <- (f(x + e) - f(lo)) / (2 * hk)
        }
        J
    }
    res <- lapply(roots, function(y) {
        p <- p_mean_given_y(y, pyr, pbr, pars)
        ev <- eigen(jac(y, p), only.values = TRUE)$values
        ev <- ev[order(-Re(ev))]
        data.frame(y = y, p = p, stable = all(Re(ev) < 0),
                   ev1_re = Re(ev[1]), ev1_im = Im(ev[1]),
                   ev2_re = Re(ev[2]), ev2_im = Im(ev[2]))
    })
    do.call(rbind, res)
}
