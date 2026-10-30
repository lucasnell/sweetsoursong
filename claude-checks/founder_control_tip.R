# R reproduction of founder_control_tip.wl: the (P_R, m_B) point where
# R_Y = 1 and R_B = 1 cross, which is the tip of the region where neither
# species can invade (regional founder control). Below the tip's m_B, the
# R_Y = 1 threshold lies above the R_B = 1 threshold, so neither invades
# between them; above it, both invade between them. Methods state that
# founder control requires m_B < 0.0091. That the region is a single wedge
# below the tip over 0.1 <= P_R <= 10 comes from the grid scan
# (founder_control_scan.wl), which r-port/check_founder_grid.R reproduces
# in R at all 7381 points. Default parameters otherwise, no Chop.
# Run from the project root (about 1 minute).
# Output: claude-checks/founder_control_tip_output.txt

suppressPackageStartupMessages(library(sweetsoursong))
options(sweetsoursong.chop = FALSE)
pars <- ss_params()

gap <- function(mb) {
    p <- update_params(pars, m_b = mb)
    inv_threshold("y", c(1.2, 1.9), p) - inv_threshold("b", c(1.2, 1.9), p)
}
mb_tip <- stats::uniroot(gap, c(0.0085, 0.0100), tol = 1e-12)$root
p_tip <- update_params(pars, m_b = mb_tip)
pr_tip <- inv_threshold("y", c(1.2, 1.9), p_tip)

out <- c(sprintf("founder_control_tip.R, %s, sweetsoursong %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M"),
                 as.character(packageVersion("sweetsoursong"))),
         sprintf("threshold gap (R_Y = 1 minus R_B = 1) at m_B = 0.0085: %.6f, at 0.0100: %.6f",
                 gap(0.0085), gap(0.0100)),
         sprintf("tip: P_R = %.10f, m_B = %.10f", pr_tip, mb_tip),
         sprintf("check: R_Y = %.10f, R_B = %.10f", inv_y(pr_tip, p_tip),
                 inv_b(pr_tip, p_tip)),
         "Mathematica (founder_control_tip_output.m): P_R = 1.5276369729, m_B = 0.0090761076")
writeLines(out, "claude-checks/founder_control_tip_output.txt")
cat(out, sep = "\n")
