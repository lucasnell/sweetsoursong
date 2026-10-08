# Fig 1E question (2026-10-08): if row 3 (P immigration) had rate
# c (m P_0R + m P_YR + m_B P_BR) / N and rows 5-6 (regional colonization)
# changed only Y or B, regional colonization would no longer add a
# pollinator in the same jump. Same mean rates; different chain. This
# compares invasion thresholds of the reduced model (vacancy closure) under
# the computed model (joint jump) and the alternative (separate events).
# In the vacancy chain the alternative adds the colonizing arrivals to the
# non-filling pollinator arrivals U_p. Run from the project root.
# Output: claude-checks/fig1e_row3_alternative_output.txt

suppressPackageStartupMessages(library(sweetsoursong))
p <- ss_params()
orig <- sweetsoursong:::vacancy_system
alt <- function(y, b, pyr, pbr, pars = ss_params()) {
    sys <- orig(y, b, pyr, pbr, pars)
    e <- pars$n - y - b
    below <- sys$p < pars$pmax
    extra <- ifelse(below, pars$c * (pars$m * pars$e_y * pyr +
                                         pars$m_b * pars$e_b * pbr) * e / pars$n^2, 0)
    sys$p_up <- sys$p_up + extra
    sys$total <- sys$total + extra
    k <- length(sys$p)
    for (i in seq_len(k)) if (i < k && sys$matrix[i, i] != 1) {
        sys$matrix[i, i] <- sys$total[i]
        sys$matrix[i, i + 1] <- -sys$p_up[i]
    }
    sys
}
thr <- function() c(y = inv_threshold("y", c(0.6, 0.85), p),
                    b = inv_threshold("b", c(2.7, 3.0), p))
t_orig <- thr()
assignInNamespace("vacancy_system", alt, "sweetsoursong")
t_alt <- thr()
assignInNamespace("vacancy_system", orig, "sweetsoursong")
lines <- c(paste0("fig1e_row3_alternative.R, ", format(Sys.time(), "%Y-%m-%d %H:%M")),
           sprintf("computed model (joint jump):  R_Y = 1 at P_R = %.6f, R_B = 1 at P_R = %.6f", t_orig[1], t_orig[2]),
           sprintf("alternative (separate events): R_Y = 1 at P_R = %.6f, R_B = 1 at P_R = %.6f", t_alt[1], t_alt[2]),
           sprintf("difference: %.2g (R_Y), %.2g (R_B)", t_alt[1] - t_orig[1], t_alt[2] - t_orig[2]))
writeLines(lines, "claude-checks/fig1e_row3_alternative_output.txt")
cat(lines, sep = "\n")
