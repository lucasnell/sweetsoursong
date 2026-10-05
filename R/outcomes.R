#' Competitive outcome from invasion criteria
#'
#' @param inv_y,inv_b Invasion criteria.
#' @return Factor with levels `"bacteria"` (only bacteria can invade),
#'   `"coexist"` (both), `"yeast"` (only yeast), and `"neither"`.
#' @export
classify_outcome <- function(inv_y, inv_b) {
    out <- ifelse(inv_y < 1 & inv_b > 1, "bacteria",
           ifelse(inv_y > 1 & inv_b > 1, "coexist",
           ifelse(inv_y > 1 & inv_b < 1, "yeast", "neither")))
    factor(out, levels = c("bacteria", "coexist", "yeast", "neither"))
}

#' Invasion criteria over a grid of regional pollinator abundance and m_b
#'
#' @param pr Values of regional pollinator abundance.
#' @param m_b Values of `m_b`.
#' @param pars An [ss_params()] object (its `m_b` is replaced).
#' @param method Passed to [inv_b()].
#' @param cores Number of cores for [parallel::mclapply()].
#' @return A data frame with `pr`, `m_b`, `inv_y`, `inv_b`, `outcome`.
#' @export
outcome_grid <- function(pr, m_b, pars = ss_params(),
                         method = c("chris", "direct"), cores = 1L) {
    method <- match.arg(method)
    grid <- expand.grid(pr = pr, m_b = m_b)
    one <- function(i) {
        p <- update_params(pars, m_b = grid$m_b[i])
        c(inv_y(grid$pr[i], p), inv_b(grid$pr[i], p, method))
    }
    res <- parallel::mclapply(seq_len(nrow(grid)), one, mc.cores = cores)
    res <- do.call(rbind, res)
    grid$inv_y <- res[, 1]
    grid$inv_b <- res[, 2]
    grid$outcome <- classify_outcome(grid$inv_y, grid$inv_b)
    grid
}
