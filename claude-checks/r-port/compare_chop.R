# How much do the manuscript results change without Chop?
# Recomputes the quantities behind Figs 3-6 with sweetsoursong.chop TRUE
# (Chris's Mathematica code) and FALSE (package default from v2.0.0).
# Run from the project root: Rscript claude-checks/r-port/compare_chop.R
# Output: claude-checks/r-port/compare_chop_output.txt

suppressPackageStartupMessages(library(sweetsoursong))
p <- ss_params()
cores <- max(1L, parallel::detectCores() - 2L)
out <- file("claude-checks/r-port/compare_chop_output.txt", "w")
say <- function(...) {
    txt <- paste0(...)
    cat(txt, "\n")
    writeLines(txt, out)
}
say("compare_chop.R, ", format(Sys.time(), "%Y-%m-%d %H:%M"), ", ",
    R.version.string, ", sweetsoursong ", as.character(packageVersion("sweetsoursong")))

run_both <- function(f) {
    options(sweetsoursong.chop = TRUE); a <- f()
    options(sweetsoursong.chop = FALSE); b <- f()
    list(chop = a, nochop = b)
}

# ---- Figs 3-4: examples ----
ex <- data.frame(pr = c(0.6, 1, 1.8, 2.6, 3), start = c(25, 25, 85, 130, 130))
r34 <- run_both(function() {
    do.call(rbind, lapply(seq_len(nrow(ex)), function(i) {
        x <- meta_equilibrium(ex$pr[i], ex$start[i], p)
        j <- joint_distribution(x, p$n * ex$pr[i] - x, p)
        d <- decompose_distribution(j)
        data.frame(pr = ex$pr[i], pyr = x, w1 = d$weights[1],
                   w2 = if (length(d$weights) > 1) d$weights[2] else NA)
    }))
})
say("\n== Figs 3-4: equilibrium P_YR and mode weights ==")
d34 <- cbind(r34$chop, nochop_pyr = r34$nochop$pyr, nochop_w1 = r34$nochop$w1)
capture <- utils::capture.output(print(d34, digits = 10))
for (l in capture) say(l)
say("max |diff P_YR| = ", signif(max(abs(r34$chop$pyr - r34$nochop$pyr)), 3),
    "; max |diff weight| = ", signif(max(abs(r34$chop$w1 - r34$nochop$w1)), 3))

# ---- Fig 5: continuation ----
cont <- function() {
    x1 <- meta_equilibrium(1, 25, p)
    run <- function(prs, x) vapply(prs, function(pr) {
        x <<- meta_equilibrium(pr, x, p)
        x
    }, 0)
    dn <- seq(0.99, 0.74, by = -0.01)
    up <- seq(1, 2.89, by = 0.01)
    v <- c(run(dn, x1), run(up, x1))
    v[order(c(dn, up))]
}
r5 <- run_both(cont)
say("\n== Fig 5: continuation, 216 values of P_R from 0.74 to 2.89 ==")
say("max |diff P_YR| = ", signif(max(abs(r5$chop - r5$nochop)), 3),
    " (max P_YR ", signif(max(r5$chop), 4), ")")

# ---- invasion thresholds (reduced model) ----
r_thr <- run_both(function() c(y = inv_threshold("y", c(0.6, 0.85), p),
                               b = inv_threshold("b", c(2.8, 3.0), p)))
say("\n== Reduced-model thresholds (eps = 1e-5) ==")
say("chop:   R_Y = 1 at ", format(r_thr$chop[["y"]], digits = 10),
    ", R_B = 1 at ", format(r_thr$chop[["b"]], digits = 10))
say("nochop: R_Y = 1 at ", format(r_thr$nochop[["y"]], digits = 10),
    ", R_B = 1 at ", format(r_thr$nochop[["b"]], digits = 10))

# ---- Fig 6: outcome grids ----
prs <- 10^seq(-1, 1, length.out = 121)
grid_fun <- function(mbs) function() outcome_grid(prs, mbs, p, cores = cores)
for (nm in c("A", "B")) {
    mbs <- if (nm == "A") seq(p$m, 0.1, length.out = 46) else seq(p$m, 1, length.out = 100)
    g <- run_both(grid_fun(mbs))
    diff <- g$chop$outcome != g$nochop$outcome
    say("\n== Fig 6", nm, ": ", nrow(g$chop), " grid points ==")
    say("outcome differs at ", sum(diff), " points")
    if (any(diff)) {
        cap <- utils::capture.output(print(cbind(g$chop[diff, c("pr", "m_b", "outcome")],
                                                 nochop = g$nochop$outcome[diff])))
        for (l in cap) say(l)
    }
    tab <- utils::capture.output(print(table(chop = g$chop$outcome,
                                             nochop = g$nochop$outcome)))
    for (l in tab) say(l)
}
close(out)
