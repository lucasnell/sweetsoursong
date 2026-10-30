# Pmax cutoff: does truncating the pollinator count at Pmax = 12 affect the
# results? For the reduced model (vacancy closure, used in the main text):
# invasion thresholds at Pmax = 12, 18, and 24; at the closed-metacommunity
# coexistence equilibria for 0.8 <= P_R <= 2.8, the stationary probability
# of P = Pmax and the relative change in the equilibrium pools (which equal
# the mean pollinator outputs E[YP] and E[BP]) from Pmax = 12 to 18. For the
# full model, the same tail probability and pool change at P_R = 1, 1.8,
# and 2.6 (Fig 3-4 examples). Equilibria are followed by continuation from
# P_R = 1 (start 25), as for Fig 5. Default parameters, no Chop. Run from
# the project root (about 10 minutes).
# Output: claude-checks/pmax_check_output.txt

suppressPackageStartupMessages(library(sweetsoursong))
options(sweetsoursong.chop = FALSE)
pars <- ss_params()
n <- pars$n
t0 <- Sys.time()

# ---- reduced model: thresholds ----
thr <- t(sapply(c(12, 18, 24), function(pm) {
    p <- update_params(pars, pmax = pm)
    c(pmax = pm, y = inv_threshold("y", c(0.6, 0.85), p),
      b = inv_threshold("b", c(2.8, 3.0), p))
}))

# ---- reduced model: equilibria along P_R ----
continuation <- function(p) {
    prs_up <- seq(1, 2.8, by = 0.1)
    prs_dn <- seq(0.9, 0.8, by = -0.1)
    run <- function(prs, x) {
        out <- numeric(length(prs))
        for (i in seq_along(prs)) {
            out[i] <- meta_equilibrium(prs[i], x, p)
            x <- out[i]
        }
        out
    }
    x1 <- meta_equilibrium(1, 25, p)
    pr <- c(rev(prs_dn), prs_up)
    pyr <- c(rev(run(prs_dn, x1)), x1, run(prs_up[-1], x1))
    data.frame(pr = pr, pyr = pyr, pbr = n * pr - pyr)
}
e12 <- continuation(pars)
e18 <- continuation(update_params(pars, pmax = 18))
stopifnot(all(e12$pyr > 1e-3 * n * e12$pr), all(e12$pbr > 1e-3 * n * e12$pr))
e12$tail_12 <- mapply(function(a, b) {
    j <- joint_distribution(a, b, pars)
    sum(j[, ncol(j)])
}, e12$pyr, e12$pbr)
e12$rel_change_pyr <- e18$pyr / e12$pyr - 1
e12$rel_change_pbr <- e18$pbr / e12$pbr - 1

# ---- full model at the Fig 3-4 coexistence examples ----
full_eq <- function(pr, start, pm) {
    p <- update_params(pars, pmax = pm)
    fe <- full_meta_equilibrium(pr, unname(start), p)
    arr <- full_distribution(n * pr - sum(fe), fe[["pyr"]], fe[["pbr"]], p)
    c(fe, tail = sum(arr[, , pm + 1]))
}
full <- do.call(rbind, lapply(c(1, 1.8, 2.6), function(pr) {
    s <- 0.99 * unlist(e12[abs(e12$pr - pr) < 1e-9, c("pyr", "pbr")])
    f12 <- full_eq(pr, s, 12)
    f18 <- full_eq(pr, f12[1:2], 18)
    data.frame(pr = pr, pyr_12 = f12[["pyr"]], pbr_12 = f12[["pbr"]],
               tail_12 = f12[["tail"]], tail_18 = f18[["tail"]],
               rel_change_pyr = f18[["pyr"]] / f12[["pyr"]] - 1,
               rel_change_pbr = f18[["pbr"]] / f12[["pbr"]] - 1)
}))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

out <- c(sprintf("pmax_check.R, %s, sweetsoursong %s, %.1f min",
                 format(Sys.time(), "%Y-%m-%d %H:%M"),
                 as.character(packageVersion("sweetsoursong")), elapsed),
         "", "== Reduced model: P_R at which R_Y = 1 and R_B = 1 ==",
         capture.output(print(as.data.frame(thr), digits = 9, row.names = FALSE)),
         "", "== Reduced model: coexistence equilibria, Pmax = 12 ==",
         "tail_12 = P(P = 12); rel_change = pool at Pmax 18 / pool at Pmax 12 - 1",
         capture.output(print(e12, digits = 4, row.names = FALSE)),
         sprintf("max tail_12 = %.3g; max |rel_change| = %.3g",
                 max(e12$tail_12),
                 max(abs(c(e12$rel_change_pyr, e12$rel_change_pbr)))),
         "", "== Full model: Fig 3-4 coexistence examples ==",
         capture.output(print(full, digits = 4, row.names = FALSE)),
         sprintf("max tail_12 = %.3g; max |rel_change| = %.3g",
                 max(full$tail_12),
                 max(abs(c(full$rel_change_pyr, full$rel_change_pbr)))))
writeLines(out, "claude-checks/pmax_check_output.txt")
cat(out, sep = "\n")
