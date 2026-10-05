#' Gillespie simulation of the reduced model with the vacancy closure
#'
#' Direct-method simulation of one plant (Chris Klausmeier's `GillespieSSA`
#' with the `proc` built by `BuildCTMCClosureTM`). Yeast-dominated flowers
#' change at the rates from [vacancy_hazards()], which depend only on Y.
#' Pollinators arrive and leave at the reduced-model rates and are not
#' truncated at `pmax`.
#'
#' Unlike Chris's code, each recorded row is the state just after a jump,
#' at the time of that jump.
#'
#' @param y0,p0 Initial state.
#' @param t_max Simulation end time.
#' @param pyr,pbr Regional pools.
#' @param pars An [ss_params()] object.
#' @param max_steps Maximum number of events.
#' @param hazards Optional output of [vacancy_hazards()] for these pools.
#' @return A data frame with columns `t`, `y`, `b`, `p`; the last row is at
#'   `t_max`.
#' @export
gillespie_reduced <- function(y0, p0, t_max, pyr, pbr, pars = ss_params(),
                              max_steps = 1e7,
                              hazards = vacancy_hazards(pyr, pbr, pars)) {
    n <- pars$n
    up <- hazards$up
    down <- hazards$down
    p_in <- pars$c * (pars$m * pyr / n + pars$m_b * pbr / n)
    y <- as.integer(y0)
    p <- as.integer(p0)
    t <- 0
    # grow output in chunks
    cap <- 10000L
    out <- matrix(NA_real_, cap, 3)
    k <- 1L
    out[k, ] <- c(t, y, p)
    steps <- 0
    repeat {
        r <- c(pars$c * p * (pars$m * y / n + pars$m_b * (n - y) / n),
               up[y + 1], down[y + 1], p_in)
        tot <- sum(r)
        if (tot <= 0 || steps >= max_steps) break
        t <- t + stats::rexp(1) / tot
        if (t >= t_max) break
        ev <- sample.int(4L, 1L, prob = r)
        switch(ev,
               p <- p - 1L,
               y <- y + 1L,
               y <- y - 1L,
               p <- p + 1L)
        steps <- steps + 1
        k <- k + 1L
        if (k > cap) {
            out <- rbind(out, matrix(NA_real_, cap, 3))
            cap <- 2L * cap
        }
        out[k, ] <- c(t, y, p)
    }
    k <- k + 1L
    if (k > cap) out <- rbind(out, matrix(NA_real_, 1, 3))
    out[k, ] <- c(t_max, y, p)
    out <- out[seq_len(k), , drop = FALSE]
    data.frame(t = out[, 1], y = out[, 2], b = n - out[, 2], p = out[, 3])
}

#' Time-weighted mean of a Gillespie path
#'
#' @param sim Output of [gillespie_reduced()].
#' @return Named vector of time-averaged `y`, `b`, `p`.
#' @export
temporal_mean <- function(sim) {
    dt <- diff(sim$t)
    k <- seq_along(dt)
    vapply(c("y", "b", "p"), function(v) sum(sim[[v]][k] * dt) / sum(dt), 0)
}
