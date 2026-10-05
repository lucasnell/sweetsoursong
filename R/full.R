# States of the full model: (y, b, p) with y + b <= n and p = 0..pmax, in
# Chris Klausmeier's order (p outermost, then b, then y).
full_states <- function(pars) {
    n <- pars$n
    yb <- do.call(rbind, lapply(0:n, function(b) cbind(y = 0:(n - b), b = b)))
    st <- do.call(rbind, lapply(0:pars$pmax, function(p) cbind(yb, p = p)))
    st
}

full_index <- function(y, b, p, n) {
    as.integer(1 + y + b * (n + 3 / 2) - b^2 / 2 + p * (n + 2) * (n + 1) / 2)
}

#' Transition-rate matrix of the full per-plant model
#'
#' Continuous-time Markov chain over (Y, B, P) with Y + B <= n and
#' P <= `pmax` (Table 1; Chris Klausmeier's `SetFullTM`). Pollinator
#' arrivals are blocked at `pmax`. The matrix is the transpose of the
#' generator: entry `[to, from]` is the rate from `from` to `to`, and columns
#' sum to zero.
#'
#' @param p0r,pyr,pbr Regional pools.
#' @param pars An [ss_params()] object.
#' @return A sparse matrix (class `dgCMatrix`) with the states as attribute
#'   `states`.
#' @export
full_tm <- function(p0r, pyr, pbr, pars = ss_params()) {
    n <- pars$n
    pmax <- pars$pmax
    st <- full_states(pars)
    y <- st[, "y"]; b <- st[, "b"]; p <- st[, "p"]
    r <- flower_rates(y, b, p, p0r, pyr, pbr, pars)
    from <- full_index(y, b, p, n)
    e <- n - y - b
    tr <- list(
        list(ok = p >= 1, dy = 0, db = 0, dp = -1, rate = "p_minus"),
        list(ok = e >= 1, dy = 1, db = 0, dp = 0, rate = "y_plus"),
        list(ok = e >= 1, dy = 0, db = 1, dp = 0, rate = "b_plus"),
        list(ok = p < pmax, dy = 0, db = 0, dp = 1, rate = "p_plus"),
        list(ok = p < pmax & e >= 1, dy = 1, db = 0, dp = 1, rate = "py_plus"),
        list(ok = p < pmax & e >= 1, dy = 0, db = 1, dp = 1, rate = "pb_plus"),
        list(ok = y >= 1, dy = -1, db = 0, dp = 0, rate = "y_minus"),
        list(ok = b >= 1, dy = 0, db = -1, dp = 0, rate = "b_minus")
    )
    ii <- jj <- integer(0)
    xx <- numeric(0)
    for (t in tr) {
        k <- which(t$ok)
        ii <- c(ii, full_index(y[k] + t$dy, b[k] + t$db, p[k] + t$dp, n))
        jj <- c(jj, from[k])
        xx <- c(xx, r[k, t$rate])
    }
    # diagonal: all eight rates below pmax; at pmax no arrivals
    at_max <- p == pmax
    diag_rate <- rowSums(r)
    diag_rate[at_max] <- rowSums(r[at_max, c("p_minus", "y_plus", "b_plus",
                                             "y_minus", "b_minus"), drop = FALSE])
    ii <- c(ii, from)
    jj <- c(jj, from)
    xx <- c(xx, -diag_rate)
    tm <- Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                               dims = rep(nrow(st), 2))
    attr(tm, "states") <- st
    tm
}

# Monoculture chain (one microbe absent) over (x, p), x = 0..n.
mono_tm <- function(which, pool, pr, pars) {
    n <- pars$n
    pmax <- pars$pmax
    st <- as.matrix(expand.grid(x = 0:n, p = 0:pmax))
    x <- st[, "x"]; p <- st[, "p"]
    p0r <- n * pr - pool
    if (which == "y") {
        r <- flower_rates(x, 0, p, p0r, pyr = pool, pbr = 0, pars)
        grow <- "y_plus"; grow_imm <- "py_plus"; die <- "y_minus"
    } else {
        r <- flower_rates(0, x, p, p0r, pyr = 0, pbr = pool, pars)
        grow <- "b_plus"; grow_imm <- "pb_plus"; die <- "b_minus"
    }
    idx <- function(x, p) as.integer(1 + x + p * (n + 1))
    from <- idx(x, p)
    tr <- list(
        list(ok = p >= 1, dx = 0, dp = -1, rate = "p_minus"),
        list(ok = x < n, dx = 1, dp = 0, rate = grow),
        list(ok = p < pmax, dx = 0, dp = 1, rate = "p_plus"),
        list(ok = p < pmax & x < n, dx = 1, dp = 1, rate = grow_imm),
        list(ok = x >= 1, dx = -1, dp = 0, rate = die)
    )
    ii <- jj <- integer(0)
    xx <- numeric(0)
    for (t in tr) {
        k <- which(t$ok)
        ii <- c(ii, idx(x[k] + t$dx, p[k] + t$dp))
        jj <- c(jj, from[k])
        xx <- c(xx, r[k, t$rate])
    }
    used <- c("p_minus", grow, "p_plus", grow_imm, die)
    diag_rate <- rowSums(r[, used])
    at_max <- p == pmax
    diag_rate[at_max] <- rowSums(r[at_max, c("p_minus", grow, die), drop = FALSE])
    ii <- c(ii, from); jj <- c(jj, from); xx <- c(xx, -diag_rate)
    tm <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = rep(nrow(st), 2))
    attr(tm, "states") <- st
    tm
}

# Solve TM v = 0 with equation j replaced by v[j] = 1. A unit row keeps the
# matrix sparse; a row of ones (sum(v) = 1) would make the LU factorization
# about 100 times slower.
tm_solve_fixed <- function(tt, j) {
    k <- tt@Dim[1]
    keep <- tt@i != (j - 1L)
    a <- Matrix::sparseMatrix(i = c(tt@i[keep] + 1L, j), j = c(tt@j[keep] + 1L, j),
                              x = c(tt@x[keep], 1), dims = c(k, k))
    f <- Matrix::lu(a, order = TRUE)
    as.numeric(Matrix::solve(f, replace(numeric(k), j, 1)))
}

# Solve TM v = 0 with sum(v) = 1 replacing the first equation. Robust but
# slow for large sparse matrices (dense row).
tm_solve_ones <- function(tm) {
    a <- tm
    a[1, ] <- 1
    f <- Matrix::lu(a, order = TRUE)
    as.numeric(Matrix::solve(f, c(1, numeric(nrow(tm) - 1))))
}

# Stationary distribution from a transposed generator (columns sum to 0).
# Fixing v[j] = 1 needs a state j that is not vanishingly rare, so the first
# pass tries states with the largest inflow / outflow ratio; the second pass
# fixes the most probable state from the first. Falls back to the dense-row
# solve if those fail. Then applies Chop the way Chris's code does: to the
# eigenvector scaled to unit Euclidean length, before normalizing to sum 1.
tm_stationary <- function(tm) {
    tt <- methods::as(tm, "TsparseMatrix")
    out_rate <- -Matrix::diag(tm)
    in_rate <- Matrix::rowSums(tm) + out_rate
    score <- in_rate / pmax(out_rate, .Machine$double.xmin)
    try_fixed <- function(j) {
        v <- tryCatch(tm_solve_fixed(tt, j), error = function(e) NULL)
        if (is.null(v) || !all(is.finite(v)) || any(v < -1e-8 * max(v))) NULL else v
    }
    v <- NULL
    for (j in utils::head(order(score, decreasing = TRUE), 5)) {
        v <- try_fixed(j)
        if (!is.null(v)) break
    }
    if (!is.null(v)) {
        j2 <- which.max(v)
        v2 <- try_fixed(j2)
        if (!is.null(v2)) v <- v2
    } else {
        v <- tm_solve_ones(tm)
    }
    v <- v / sum(v)
    if (use_chop()) {
        u <- v / sqrt(sum(v^2))
        u <- chop(u)
        v <- u / sum(u)
    } else {
        # negative entries are round-off of order 1e-17
        v <- pmax(v, 0)
        v <- v / sum(v)
    }
    v
}

#' Stationary distribution of the full per-plant model
#'
#' @inheritParams full_tm
#' @return A 3-D array of probabilities indexed by y, b, p (zero where
#'   y + b > n).
#' @export
full_distribution <- function(p0r, pyr, pbr, pars = ss_params()) {
    tm <- full_tm(p0r, pyr, pbr, pars)
    v <- tm_stationary(tm)
    st <- attr(tm, "states")
    n <- pars$n
    arr <- array(0, c(n + 1, n + 1, pars$pmax + 1),
                 dimnames = list(y = 0:n, b = 0:n, p = 0:pars$pmax))
    arr[cbind(st[, "y"] + 1, st[, "b"] + 1, st[, "p"] + 1)] <- v
    arr
}

#' Pollinator outputs of the full model
#'
#' `full_in_out()` gives `E[Y P]` and `E[B P]` for given pools, with
#' `p0r = n * pr - pyr - pbr`. `full_mono_out()` gives `E[X P]` for a
#' monoculture of species `which` with pool `pool`.
#'
#' @param pyr,pbr Regional pools.
#' @param pr Regional pollinator abundance.
#' @param pars An [ss_params()] object.
#' @export
full_in_out <- function(pyr, pbr, pr, pars = ss_params()) {
    arr <- full_distribution(pars$n * pr - pyr - pbr, pyr, pbr, pars)
    g <- full_grids(arr)
    c(yp = sum(arr * g$y * g$p), bp = sum(arr * g$b * g$p))
}

full_grids <- function(arr) {
    d <- dim(arr)
    list(y = slice.index(arr, 1) - 1, b = slice.index(arr, 2) - 1,
         p = slice.index(arr, 3) - 1)
}

#' @rdname full_in_out
#' @param which `"y"` or `"b"`.
#' @param pool The monoculture species' regional pool.
#' @export
full_mono_out <- function(which = c("y", "b"), pool, pr, pars = ss_params()) {
    which <- match.arg(which)
    tm <- mono_tm(which, pool, pr, pars)
    v <- tm_stationary(tm)
    st <- attr(tm, "states")
    sum(v * st[, "x"] * st[, "p"])
}

#' Invasion criteria of the full model
#'
#' The resident's monoculture equilibrium solves
#' `full_mono_out(resident, pool) = pool` from `0.99 * n * pr`; the invader's
#' criterion is its pollinator output divided by `eps` (Chris Klausmeier's
#' `FullInvY`, `FullInvB`).
#'
#' @inheritParams full_in_out
#' @export
full_inv_y <- function(pr, pars = ss_params()) {
    pbr <- newton_root(function(x) full_mono_out("b", x, pr, pars) - x,
                       0.99 * pars$n * pr, lower = 0, upper = pars$n * pr)
    full_in_out(pars$eps, pbr, pr, pars)[["yp"]] / pars$eps
}

#' @rdname full_inv_y
#' @export
full_inv_b <- function(pr, pars = ss_params()) {
    pyr <- newton_root(function(x) full_mono_out("y", x, pr, pars) - x,
                       0.99 * pars$n * pr, lower = 0, upper = pars$n * pr)
    full_in_out(pyr, pars$eps, pr, pars)[["bp"]] / pars$eps
}

#' Closed-metacommunity equilibrium of the full model
#'
#' Solves `full_in_out(pyr, pbr) = (pyr, pbr)` by Newton's method in two
#' dimensions from `start`.
#'
#' @inheritParams full_in_out
#' @param start Starting values `c(pyr, pbr)`.
#' @param tol Convergence tolerance.
#' @return Named vector `c(pyr, pbr)`.
#' @export
full_meta_equilibrium <- function(pr, start, pars = ss_params(), tol = 1e-10) {
    f <- function(x) full_in_out(x[1], x[2], pr, pars) - x
    x <- start
    fx <- f(x)
    for (it in 1:100) {
        J <- matrix(0, 2, 2)
        for (k in 1:2) {
            h <- 1e-6 * max(1, abs(x[k]))
            e <- c(0, 0); e[k] <- h
            J[, k] <- (f(x + e) - f(x - e)) / (2 * h)
        }
        step <- -solve(J, fx)
        lam <- 1
        repeat {
            xn <- x + lam * step
            if (all(xn > 0) && sum(xn) < pars$n * pr) {
                fn <- f(xn)
                if (sum(fn^2) < sum(fx^2) || lam < 1e-8) break
            }
            lam <- lam / 2
            if (lam < 1e-12) stop("Newton step failed")
        }
        if (max(abs(xn - x)) < tol * max(1, abs(x))) {
            x <- xn
            break
        }
        x <- xn; fx <- fn
    }
    c(pyr = x[1], pbr = x[2])
}
