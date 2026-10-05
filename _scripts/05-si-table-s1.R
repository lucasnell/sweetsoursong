# SI Table S1: values of P_R at which each invasion criterion equals 1, for
# the reduced model (vacancy closure) and the full model, across epsilon and
# Pmax. Table S1 was computed without Chop (Chop sets InvY to 0 in the
# reduced model when epsilon <= 1e-7), so this script turns it off.

source("_scripts/00-preamble.R")
options(sweetsoursong.chop = FALSE)

thr_reduced <- function(eps) {
    p <- update_params(pars, eps = eps)
    c(y = inv_threshold("y", c(0.6, 0.85), p),
      b = inv_threshold("b", c(2.8, 3.0), p, method = "direct"))
}
thr_full <- function(eps, pmax) {
    p <- update_params(pars, eps = eps, pmax = pmax)
    c(y = stats::uniroot(function(pr) full_inv_y(pr, p) - 1, c(0.7, 0.76),
                         tol = 1e-9)$root,
      b = stats::uniroot(function(pr) full_inv_b(pr, p) - 1, c(2.8, 2.9),
                         tol = 1e-9)$root)
}

tab <- cached("si-table-s1.rds", {
    rows <- list(
        data.frame(model = "reduced", eps = 1e-5, pmax = 12, t(thr_reduced(1e-5))),
        data.frame(model = "reduced", eps = 1e-7, pmax = 12, t(thr_reduced(1e-7))),
        data.frame(model = "reduced", eps = 1e-9, pmax = 12, t(thr_reduced(1e-9))),
        data.frame(model = "full", eps = 1e-5, pmax = 12, t(thr_full(1e-5, 12))),
        data.frame(model = "full", eps = 1e-7, pmax = 12, t(thr_full(1e-7, 12))),
        data.frame(model = "full", eps = 1e-5, pmax = 18, t(thr_full(1e-5, 18))),
        data.frame(model = "full", eps = 1e-7, pmax = 18, t(thr_full(1e-7, 18))))
    do.call(rbind, rows)
})
names(tab)[4:5] <- c("pr_inv_y_eq_1", "pr_inv_b_eq_1")
write.csv(tab, data_path("si-table-s1.csv"), row.names = FALSE)
print(tab, digits = 7)
