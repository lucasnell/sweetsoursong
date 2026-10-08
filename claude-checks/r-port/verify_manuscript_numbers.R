# Compare the R port's outputs with numbers printed in the manuscript
# (sweetsoursong-ms, origin/main on 2026-10-05). Run from the project root
# after the _scripts/ figure scripts.
# Inputs: _data/si-table-s2.csv, _data/fig3-4-equilibria.csv,
#         _data/fig5-continuation.csv, _data/fig6-grid-*.rds
# Output: claude-checks/r-port/verify_manuscript_numbers_output.txt

suppressPackageStartupMessages(library(sweetsoursong))
lines <- paste0("verify_manuscript_numbers.R, ", format(Sys.time(), "%Y-%m-%d %H:%M"),
                ", sweetsoursong ", as.character(packageVersion("sweetsoursong")))
add <- function(...) lines <<- c(lines, paste0(...))
check <- function(what, r_value, ms_value, digits) {
    ok <- round(r_value, digits) == round(ms_value, digits)
    add(sprintf("%-55s R %-12s MS %-10s %s", what, format(r_value, digits = 10),
                format(ms_value), if (ok) "match" else "DIFFERENT"))
}

# SI Table S2 (99-supplement.tex, tab:numerics), digits as printed
s1 <- read.csv("_data/si-table-s2.csv")
add("\n== SI Table S2 ==")
ms <- data.frame(model = c("reduced", "reduced", "reduced", "full", "full"),
                 eps = c(1e-5, 1e-7, 1e-9, 1e-5, 1e-5), pmax = c(12, 12, 12, 12, 18),
                 y = c(0.730855, 0.730855, 0.730855, 0.729677, 0.729670),
                 b = c(2.895810, 2.895811, 2.89581, 2.867489, 2.867473),
                 dig_b = c(6, 6, 5, 6, 6))
for (i in seq_len(nrow(ms))) {
    r <- s1[s1$model == ms$model[i] & s1$eps == ms$eps[i] & s1$pmax == ms$pmax[i], ]
    lab <- sprintf("%s, eps %g, Pmax %d", ms$model[i], ms$eps[i], ms$pmax[i])
    check(paste(lab, "R_Y = 1"), r$pr_inv_y_eq_1, ms$y[i], 6)
    check(paste(lab, "R_B = 1"), r$pr_inv_b_eq_1, ms$b[i], ms$dig_b[i])
}
r <- s1[s1$model == "full" & s1$eps == 1e-7, ]
add("full model at eps 1e-7 (Table S2 row 'eps 1e-5 to 1e-7'): ",
    paste(format(unlist(r[, 4:5]), digits = 7), collapse = ", "))

# Fig 3 caption (98-figures.tex): P_R = 1.8, P_YR = 71.87, P_BR = 18.13
e <- read.csv("_data/fig3-4-equilibria.csv")
add("\n== Fig 3 caption ==")
check("P_YR at P_R = 1.8", e$pyr[e$pr == 1.8], 71.87, 2)
check("P_BR at P_R = 1.8", e$pbr[e$pr == 1.8], 18.13, 2)

# Fig 4 caption: outcomes at P_R = 0.6, 1, 1.8, 2.6, 3
add("\n== Fig 4 panels ==")
for (i in seq_len(nrow(e))) {
    add(sprintf("P_R = %.1f: R_Y %s 1, R_B %s 1, modes %d, weights %s", e$pr[i],
                if (e$inv_y[i] > 1) ">" else "<", if (e$inv_b[i] > 1) ">" else "<",
                e$n_modes[i], e$weights[i]))
}

# Fig 5 caption: coexistence for 0.73 < P_R < 2.90
add("\n== Fig 5 caption (coexistence for 0.73 < P_R < 2.90) ==")
th_y <- inv_threshold("y", c(0.6, 0.85))
th_b <- inv_threshold("b", c(2.8, 3.0))
check("R_Y = 1 at P_R", th_y, 0.73, 2)
check("R_B = 1 at P_R", th_b, 2.90, 2)
f5 <- read.csv("_data/fig5-continuation.csv")
add("continuation: two modes for P_R in ", paste(range(f5$pr[f5$n_modes == 2]), collapse = " to "),
    "; R_Y > 1 and R_B > 1 at all ", nrow(f5), " points: ",
    all(f5$inv_y > 1 & f5$inv_b > 1))

# Fig 6: outcome counts; no point where neither species can invade
add("\n== Fig 6 ==")
for (nm in c("a", "b")) {
    g <- readRDS(sprintf("_data/fig6-grid-%s.rds", nm))
    tab <- table(g$outcome)
    add(sprintf("panel %s: %s", toupper(nm),
                paste(names(tab), tab, sep = " ", collapse = ", ")))
}

# ms_numbers_output.m: m_B = m window 1.46 < P_R < 1.57 (abstract, CLAUDE.md)
add("\n== No pollinator preference (m_B = m) ==")
pm <- update_params(ss_params(), m_b = 0.01)
check("R_Y = 1 at P_R", inv_threshold("y", c(1.4, 1.5), pm), 1.462006020668476, 8)
check("R_B = 1 at P_R", inv_threshold("b", c(1.5, 1.6), pm), 1.566833665276046, 8)

writeLines(lines, "claude-checks/r-port/verify_manuscript_numbers_output.txt")
cat(lines, sep = "\n")
