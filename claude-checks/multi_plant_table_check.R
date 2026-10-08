# Check that the K-plant transitions (manuscript SI Table S1; Table 2 when
# written) reproduce the one-plant transitions with a regional pool
# (manuscript Fig 1E; Table 1 when written):
#  (1) summed over destinations, departures from plant j equal the one-plant
#      emigration rate;
#  (2) summed over source plants, arrivals at plant i (by type) equal the
#      one-plant immigration rates with P_XR = sum_{j != i} P_j X_j / (K - 1).
# Run from the project root. Output: claude-checks/multi_plant_table_check_output.txt

suppressPackageStartupMessages(library(sweetsoursong))
p <- update_params(ss_params(), m_b = 0.05)
n <- p$n; c <- p$c; m <- p$m; mb <- p$m_b; ey <- p$e_y; eb <- p$e_b
set.seed(3)
lines <- c(paste0("multi_plant_table_check.R, ", format(Sys.time(), "%Y-%m-%d %H:%M")))
worst <- 0
for (rep in 1:200) {
    K <- sample(2:20, 1)
    Y <- sample(0:n, K, replace = TRUE)
    B <- vapply(Y, function(y) sample(0:(n - y), 1), 0)
    E <- n - Y - B
    P <- sample(0:8, K, replace = TRUE)
    i <- sample.int(K, 1)
    js <- setdiff(seq_len(K), i)
    # Table 2 rates for moves j -> i
    move_only <- c * P[js] * (m * E[js] + m * (1 - ey * E[i] / n) * Y[js] +
                              mb * (1 - eb * E[i] / n) * B[js]) / (n * (K - 1))
    move_y <- c * m * ey * P[js] * Y[js] * E[i] / (n^2 * (K - 1))
    move_b <- c * mb * eb * P[js] * B[js] * E[i] / (n^2 * (K - 1))
    # Table 1 at plant i with pools from the other plants
    pools <- c(p0r = sum(P[js] * E[js]), pyr = sum(P[js] * Y[js]),
               pbr = sum(P[js] * B[js])) / (K - 1)
    t1 <- flower_rates(Y[i], B[i], P[i], pools[["p0r"]], pools[["pyr"]], pools[["pbr"]], p)
    d_arr <- c(sum(move_only) - t1[, "p_plus"], sum(move_y) - t1[, "py_plus"],
               sum(move_b) - t1[, "pb_plus"])
    # departures from plant i summed over the K - 1 destinations
    dep <- c * P[i] * (m * E[i] + m * Y[i] + mb * B[i]) / n
    others <- setdiff(seq_len(K), i)
    dep2 <- sum(vapply(others, function(k) {
        c * P[i] * (m * E[i] + m * (1 - ey * E[k] / n) * Y[i] +
                    mb * (1 - eb * E[k] / n) * B[i]) / (n * (K - 1)) +
            c * m * ey * P[i] * Y[i] * E[k] / (n^2 * (K - 1)) +
            c * mb * eb * P[i] * B[i] * E[k] / (n^2 * (K - 1))
    }, 0))
    d_dep <- c(dep2 - dep, t1[, "p_minus"] - dep)
    worst <- max(worst, abs(c(d_arr, d_dep)) / pmax(1, abs(c(t1[, c("p_plus", "py_plus", "pb_plus")], dep, dep))))
}
lines <- c(lines, "200 random states, K = 2..20 plants",
           paste0("max relative difference, Table 2 summed vs Table 1: ", signif(worst, 3)))
writeLines(lines, "claude-checks/multi_plant_table_check_output.txt")
cat(lines, sep = "\n")
