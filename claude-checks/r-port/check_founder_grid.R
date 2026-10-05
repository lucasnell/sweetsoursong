# R port against the Mathematica founder-control scan: invasion criteria
# and outcomes at all 7381 (PR, mB) points of founder_control_scan.wl
# (Chop off, Chris Klausmeier's InvB). Run from the project root.
# Input: claude-checks/r-port/founder_control_grid.csv (export_founder_grid.wl)
# Output: claude-checks/r-port/check_founder_grid_output.txt

suppressPackageStartupMessages(library(sweetsoursong))
options(sweetsoursong.chop = FALSE)
p <- ss_params()
ref <- read.csv("claude-checks/r-port/founder_control_grid.csv")
cores <- max(1L, parallel::detectCores() - 2L)
res <- parallel::mclapply(seq_len(nrow(ref)), function(i) {
    pp <- update_params(p, m_b = ref$mB[i])
    c(inv_y(ref$pr[i], pp), inv_b(ref$pr[i], pp, "chris"), inv_b(ref$pr[i], pp, "direct"))
}, mc.cores = cores)
res <- do.call(rbind, res)
out_m <- classify_outcome(ref$invY, ref$invB)
out_r <- classify_outcome(res[, 1], res[, 3])
rel <- function(a, b) abs(a - b) / pmax(1, abs(b))
lines <- c(
    paste0("check_founder_grid.R, ", format(Sys.time(), "%Y-%m-%d %H:%M"),
           ", sweetsoursong ", as.character(packageVersion("sweetsoursong"))),
    paste0("points: ", nrow(ref)),
    paste0("max relative difference, R_Y: ", signif(max(rel(res[, 1], ref$invY)), 3)),
    paste0("max relative difference, R_B (Chris's form): ", signif(max(rel(res[, 2], ref$invB)), 3)),
    paste0("outcome differs (R direct R_B vs Mathematica): ", sum(out_m != out_r)),
    utils::capture.output(print(table(mathematica = out_m, r = out_r))))
writeLines(lines, "claude-checks/r-port/check_founder_grid_output.txt")
cat(lines, sep = "\n")
