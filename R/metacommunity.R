#' Expected yeast-carrying pollinator output of one plant
#'
#' `E[Y P]` over the stationary distribution of a plant given the regional
#' pools, with `pbr = n * pr - pyr` (Chris Klausmeier's `InOutPY`). In a
#' closed metacommunity, `pyr = e_yp(pyr)`.
#'
#' @param pyr Regional pool of pollinators that last visited a
#'   yeast-dominated flower.
#' @param pr Regional pollinator abundance (mean pollinators per plant).
#' @param pars An [ss_params()] object.
#' @export
e_yp <- function(pyr, pr, pars = ss_params()) {
    pbr <- pars$n * pr - pyr
    y <- 0:pars$n
    pi_y <- y_distribution(pyr, pbr, pars)
    sum(pi_y * y * p_mean_given_y(y, pyr, pbr, pars))
}

#' Expected bacteria-carrying pollinator output of one plant
#'
#' `E[B P]` with `B = n - Y`.
#'
#' @inheritParams e_yp
#' @param pbr Regional bacteria pool. Pass it directly when it is small:
#'   `n * pr - pyr` loses relative precision when `pyr` is close to `n * pr`.
#' @export
e_bp <- function(pyr, pr, pars = ss_params(), pbr = pars$n * pr - pyr) {
    y <- 0:pars$n
    pi_y <- y_distribution(pyr, pbr, pars)
    sum(pi_y * (pars$n - y) * p_mean_given_y(y, pyr, pbr, pars))
}

#' Invasion criteria in the closed metacommunity
#'
#' `inv_y()` is `E[Y P] / eps` when the yeast pool is `eps` in a
#' bacteria-dominated metacommunity. `inv_b()` is the analogue for bacteria.
#' Each species can invade when its criterion exceeds 1.
#'
#' `inv_b(method = "chris")` is `(n * pr - E[Y P]) / eps`, as in Chris
#' Klausmeier's `InvB`. It has the same threshold as `method = "direct"`,
#' `E[B P] / eps`, but a different magnitude:
#' `direct - 1 = (m / m_b) * (chris - 1)`.
#'
#' @inheritParams e_yp
#' @param method For `inv_b()`, `"chris"` or `"direct"`.
#' @export
inv_y <- function(pr, pars = ss_params()) {
    e_yp(pars$eps, pr, pars) / pars$eps
}

#' @rdname inv_y
#' @export
inv_b <- function(pr, pars = ss_params(), method = c("chris", "direct")) {
    method <- match.arg(method)
    pyr <- pars$n * pr - pars$eps
    switch(method,
           chris = (pars$n * pr - e_yp(pyr, pr, pars)) / pars$eps,
           direct = e_bp(pyr, pr, pars, pbr = pars$eps) / pars$eps)
}

#' Closed-metacommunity equilibrium
#'
#' Solves `e_yp(pyr, pr) = pyr` for `pyr` in `[0, n * pr]`. The function
#' `e_yp(pyr) - pyr` tends to 0 at both ends (one species absent).
#'
#' First, Newton's method from `start`, as Mathematica's `FindRoot` does in
#' Chris Klausmeier's notebook; its result is used if it is an interior
#' root. Otherwise (Newton reaches an end, or `start` is at an end, where
#' `e_yp` jumps because the end state has zero probability), the result is
#' the candidate nearest to `start` among the sign changes found by
#' evaluating at geometrically spaced points on both sides of `start` and
#' the ends: 0 if the function is negative just above it (yeast absent),
#' `n * pr` if positive just below it (bacteria absent).
#'
#' This reproduces the roots that `FindRoot` finds from the starting values
#' in the notebook (Figs 3-5).
#'
#' @inheritParams e_yp
#' @param start Starting value of `pyr`.
#' @param tol Tolerance on `pyr`, relative to `n * pr`.
#' @return The equilibrium `pyr`.
#' @export
meta_equilibrium <- function(pr, start, pars = ss_params(), tol = 1e-12) {
    s <- pars$n * pr
    f <- function(x) e_yp(x, pr, pars) - x
    edge <- 1e-6 * s
    x0 <- min(max(start, 1e-9 * s), s - 1e-9 * s)
    x <- tryCatch(suppressWarnings(newton_root(f, x0, 0, s, tol = tol)),
                  error = function(e) NA_real_)
    if (!is.na(x) && x > edge && x < s - edge && abs(f(x)) < 1e-8 * s) return(x)
    nearest_sign_change(f, start, 0, s, tol = tol * s)
}

# Root of f nearest to x0 in [lo, hi], where f tends to 0 at both ends.
# Evaluates f at x0 and at x0 -/+ d for d = 1e-4, 2e-4, ... times (hi - lo)
# inside the open interval. Candidates are the sign changes (refined with
# uniroot) and the ends: lo if f < 0 just above it, hi if f > 0 just below
# it (f moves toward 0 there, as at an attracting root). Returns the
# candidate nearest to x0.
nearest_sign_change <- function(f, x0, lo, hi, tol) {
    w <- hi - lo
    inset <- 1e-9 * w
    x0c <- min(max(x0, lo + inset), hi - inset)
    d <- w * 1e-4 * 2^(0:14)
    xs <- sort(unique(c(x0c, pmax(x0c - d, lo + inset), pmin(x0c + d, hi - inset))))
    fs <- vapply(xs, f, 0)
    nx <- length(xs)
    cand <- numeric(0)
    cand_dist <- numeric(0)
    zero <- xs[fs == 0]
    cand <- c(cand, zero)
    cand_dist <- c(cand_dist, abs(zero - x0))
    sc <- which(fs[-1] != 0 & fs[-nx] != 0 & sign(fs[-1]) != sign(fs[-nx]))
    cand_dist <- c(cand_dist, pmin(abs(xs[sc] - x0), abs(xs[sc + 1] - x0)))
    cand <- c(cand, -sc)  # negative: bracket index, refined below if chosen
    if (fs[1] < 0) {
        cand <- c(cand, NA)
        cand_dist <- c(cand_dist, abs(x0 - lo))
    }
    if (fs[nx] > 0) {
        cand <- c(cand, Inf)
        cand_dist <- c(cand_dist, abs(hi - x0))
    }
    if (length(cand) == 0) stop("no root found")
    k <- which.min(cand_dist)
    best <- cand[k]
    if (is.na(best)) return(lo)
    if (is.infinite(best)) return(hi)
    if (k <= length(zero)) return(best)
    i <- -best
    stats::uniroot(f, xs[c(i, i + 1)], f.lower = fs[i], f.upper = fs[i + 1],
                   tol = tol)$root
}

# Newton's method with a finite-difference derivative (one-sided at the
# bounds) and step halving to reduce |f|. Steps are clamped to
# [lower, upper], so roots on a bound (one species absent) can be reached.
newton_root <- function(f, x0, lower = -Inf, upper = Inf, tol = 1e-10,
                        maxit = 100) {
    x <- x0
    fx <- f(x)
    if (fx == 0) return(x)
    for (i in seq_len(maxit)) {
        h <- 1e-6 * max(1, abs(x))
        if (x - h < lower) {
            dfx <- (f(x + h) - fx) / h
        } else if (x + h > upper) {
            dfx <- (fx - f(x - h)) / h
        } else {
            dfx <- (f(x + h) - f(x - h)) / (2 * h)
        }
        if (!is.finite(dfx) || dfx == 0) stop("zero or non-finite derivative at ", x)
        step <- -fx / dfx
        lam <- 1
        repeat {
            xn <- min(max(x + lam * step, lower), upper)
            fn <- f(xn)
            if (abs(fn) < abs(fx) || lam < 1e-8) break
            lam <- lam / 2
        }
        if (fn == 0 || abs(xn - x) < tol * max(1, abs(x))) return(xn)
        x <- xn
        fx <- fn
    }
    warning("Newton's method did not converge")
    x
}

#' Regional pollinator abundance at which an invasion criterion equals 1
#'
#' @param which `"y"` or `"b"`.
#' @param interval Interval of `pr` containing exactly one crossing.
#' @inheritParams inv_y
#' @param tol Tolerance passed to [stats::uniroot()].
#' @export
inv_threshold <- function(which = c("y", "b"), interval, pars = ss_params(),
                          method = c("chris", "direct"), tol = 1e-10) {
    which <- match.arg(which)
    method <- match.arg(method)
    f <- switch(which,
                y = function(pr) inv_y(pr, pars) - 1,
                b = function(pr) inv_b(pr, pars, method) - 1)
    stats::uniroot(f, interval, tol = tol)$root
}
