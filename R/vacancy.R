#' Pollinator chain for one uncolonized flower
#'
#' After a flower senesces, the plant has `y` yeast-dominated flowers, `b`
#' bacteria-dominated flowers, and one uncolonized flower. Pollinators arrive
#' and leave (P = 0, ..., `pmax`) until the uncolonized flower is filled.
#' This builds that chain (Chris Klausmeier's
#' `BuildCTMCClosureVacancySystem`).
#'
#' @param y,b Numbers of yeast- and bacteria-dominated flowers after the
#'   senescence (`y + b <= n`).
#' @param pyr,pbr Regional pools.
#' @param pars An [ss_params()] object.
#' @return A list with the matrix and rate vectors over P = 0, ..., `pmax`.
#' @export
vacancy_system <- function(y, b, pyr, pbr, pars = ss_params()) {
    if (y + b > pars$n) stop("y + b must be <= n")
    pmax <- pars$pmax
    ps <- 0:pmax
    r <- flower_rates(y, b, ps, p0r = 0, pyr = pyr, pbr = pbr, pars = pars,
                      e = pars$n - y - b)
    below <- ps < pmax
    y_repl <- r[, "y_plus"] + ifelse(below, r[, "py_plus"], 0)
    b_repl <- r[, "b_plus"] + ifelse(below, r[, "pb_plus"], 0)
    p_down <- ifelse(ps > 0, r[, "p_minus"], 0)
    p_up <- ifelse(below, r[, "p_plus"], 0)
    total <- y_repl + b_repl + p_down + p_up

    k <- pmax + 1L
    mat <- matrix(0, k, k)
    for (i in seq_len(k)) {
        if (maybe_chop(total[i]) == 0) {
            mat[i, i] <- 1
        } else {
            mat[i, i] <- total[i]
            if (i > 1) mat[i, i - 1] <- -p_down[i]
            if (i < k) mat[i, i + 1] <- -p_up[i]
        }
    }
    list(p = ps, matrix = mat, y_repl = y_repl, b_repl = b_repl,
         p_down = p_down, p_up = p_up, total = total)
}

#' Probabilities that an uncolonized flower is filled by yeast or bacteria
#'
#' Absorption probabilities of the chain in [vacancy_system()], from each
#' starting pollinator count (Chris Klausmeier's
#' `CTMCClosureReplacementProbabilityVectors`).
#'
#' @inheritParams vacancy_system
#' @return A matrix with columns `h_y` and `h_b`, one row per P.
#' @export
fill_probs <- function(y, b, pyr, pbr, pars = ss_params()) {
    sys <- vacancy_system(y, b, pyr, pbr, pars)
    sol <- solve(sys$matrix, cbind(sys$y_repl, sys$b_repl))
    sol <- pmin(pmax(maybe_chop(sol), 0), 1)
    colnames(sol) <- c("h_y", "h_b")
    sol
}

#' Poisson weights for the pollinator count at senescence
#'
#' @param y Number of yeast-dominated flowers before the senescence.
#' @inheritParams vacancy_system
#' @param normalize Renormalize the Poisson probabilities truncated at
#'   `pmax` to sum to 1.
#' @export
poisson_weights <- function(y, pyr, pbr, pars = ss_params(), normalize = TRUE) {
    lambda <- p_mean_given_y(y, pyr, pbr, pars)
    w <- stats::dpois(0:pars$pmax, lambda)
    mass <- sum(w)
    if (normalize) {
        # in log space: for large lambda the truncated probabilities underflow
        # to 0 in double precision (Mathematica switches to extended
        # precision there instead)
        lw <- stats::dpois(0:pars$pmax, lambda, log = TRUE)
        w <- exp(lw - max(lw))
        w <- w / sum(w)
    }
    attr(w, "lambda") <- lambda
    attr(w, "tail_mass") <- max(0, 1 - mass)
    w
}

#' Rates of change in yeast-dominated flowers under the vacancy closure
#'
#' For each plant state y = 0, ..., n, the rate at which a bacteria-dominated
#' flower senesces and is replaced by yeast (`up`) and the reverse (`down`),
#' averaging over the pollinator count with `P | Y ~ Poisson(E[P | Y])` truncated
#' at `pmax` (Chris Klausmeier's `CTMCClosureHazards`).
#'
#' @param pyr,pbr Regional pools.
#' @param pars An [ss_params()] object.
#' @return A data frame with columns `y`, `up`, `down`, `lambda`, `tail_mass`.
#' @export
vacancy_hazards <- function(pyr, pbr, pars = ss_params()) {
    n <- pars$n
    # The vacancy after a bacteria death in state y is (y, n - y - 1); the
    # vacancy after a yeast death in state y + 1 is the same state.
    h <- lapply(0:(n - 1), function(k) fill_probs(k, n - 1 - k, pyr, pbr, pars))
    up <- down <- lambda <- tail_mass <- numeric(n + 1)
    for (y in 0:n) {
        w <- poisson_weights(y, pyr, pbr, pars)
        lambda[y + 1] <- attr(w, "lambda")
        tail_mass[y + 1] <- attr(w, "tail_mass")
        b <- n - y
        if (b > 0) up[y + 1] <- pars$d * b * sum(w * h[[y + 1]][, "h_y"])
        if (y > 0) down[y + 1] <- pars$d * y * sum(w * h[[y]][, "h_b"])
    }
    data.frame(y = 0:n, up = up, down = down, lambda = lambda,
               tail_mass = tail_mass)
}
