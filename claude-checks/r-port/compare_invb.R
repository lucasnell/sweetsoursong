# Does the Fig 6 outcome map change with the direct R_B (E[BP] / eps,
# package default from v2.0.0) instead of Chris Klausmeier's
# (n PR - E[YP]) / eps? Run from the project root.
# Output: claude-checks/r-port/compare_invb_output.txt

suppressPackageStartupMessages(library(sweetsoursong))
options(sweetsoursong.chop = FALSE)
p <- ss_params()
cores <- max(1L, parallel::detectCores() - 2L)
out <- file("claude-checks/r-port/compare_invb_output.txt", "w")
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); writeLines(txt, out) }
say("compare_invb.R, ", format(Sys.time(), "%Y-%m-%d %H:%M"), ", sweetsoursong ",
    as.character(packageVersion("sweetsoursong")))
prs <- 10^seq(-1, 1, length.out = 121)
for (nm in c("A", "B")) {
    mbs <- if (nm == "A") seq(p$m, 0.1, length.out = 46) else seq(p$m, 1, length.out = 100)
    gd <- outcome_grid(prs, mbs, p, method = "direct", cores = cores)
    gc <- outcome_grid(prs, mbs, p, method = "chris", cores = cores)
    diff <- gd$outcome != gc$outcome
    say("Fig 6", nm, ": ", nrow(gd), " points; outcome differs at ", sum(diff))
    # identity check: direct - 1 = (m / m_b) (chris - 1)
    pred <- 1 + (p$m / gd$m_b) * (gc$inv_b - 1)
    say("  max relative deviation from the identity: ",
        signif(max(abs(gd$inv_b - pred) / pmax(1, abs(gd$inv_b))), 3))
    if (any(diff)) {
        cap <- utils::capture.output(print(cbind(gd[diff, c("pr", "m_b", "inv_b")],
                                                 chris = gc$inv_b[diff])))
        for (l in cap) say(l)
    }
}
close(out)
