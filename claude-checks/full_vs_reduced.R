# Does the reduced model (no uncolonized flowers, vacancy closure) give the
# same closed-metacommunity distributions as the full model with uncolonized
# flowers? SI Table S2 compares only invasion thresholds. This compares, at
# P_R values across the coexistence range (Figs 3-5): the regional pools,
# mean Y and P, the weight of the yeast-dominated mode, and the mean Y and P
# within each mode. Modes come from the same watershed decomposition
# (decompose_distribution) applied to the (Y, P) distribution of each model;
# for the full model that is the joint distribution summed over B.
# Distances between the two distributions:
#  tv_yp: total variation distance, sum |reduced - full| / 2, between the
#    (Y, P) distributions. It is large whenever the full model moves mass
#    by one flower: with 0.3-0.6 uncolonized flowers per plant, a plant at
#    Y = 50 in the reduced model is often at Y = 49 in the full model.
#  w1_y, w1_b: Wasserstein-1 distance between the marginal distributions of
#    Y and of B (B = N - Y in the reduced model), the sum over k of
#    |F_reduced(k) - F_full(k)| for cumulative distributions F. It is the
#    mean number of flowers by which mass must be moved to turn one
#    distribution into the other.
# Default parameters (ss_params()), no Chop. Run from the project root
# (about 3 minutes).
# Output: claude-checks/full_vs_reduced_output.txt, full_vs_reduced_output.csv

suppressPackageStartupMessages(library(sweetsoursong))
options(sweetsoursong.chop = FALSE)
pars <- ss_params()
n <- pars$n

prs <- c(0.8, 1, 1.4, 1.8, 2.2, 2.6, 2.8)

summarize <- function(joint) {
    d <- decompose_distribution(joint)
    iy <- which.max(d$means[, "y"])
    ib <- which.min(d$means[, "y"])
    m <- joint_means(joint)
    data.frame(mean_y = m[["y"]], mean_p = m[["p"]],
               n_modes = length(d$weights),
               w_y = d$weights[iy],
               y_on_y = d$means[iy, "y"], p_on_y = d$means[iy, "p"],
               y_on_b = d$means[ib, "y"], p_on_b = d$means[ib, "p"])
}

one <- function(pr) {
    # reduced model; start value from the Fig 5 continuation (pyr / (n pr)
    # rises from about 0.1 to 0.9 across the coexistence range)
    start <- n * pr * min(max((pr - 0.7) / 2.3, 0.05), 0.95)
    pyr <- meta_equilibrium(pr, start, pars)
    pbr <- n * pr - pyr
    stopifnot(pyr > 1e-3 * n * pr, pbr > 1e-3 * n * pr)  # interior root
    jr <- joint_distribution(pyr, pbr, pars)
    # full model, starting from the reduced pools scaled to leave room for
    # the uncolonized-flower pool P_0R
    fe <- full_meta_equilibrium(pr, 0.99 * c(pyr, pbr), pars)
    p0r <- n * pr - sum(fe)
    arr <- full_distribution(p0r, fe[["pyr"]], fe[["pbr"]], pars)
    jf <- apply(arr, c(1, 3), sum)
    dimnames(jf) <- list(y = 0:n, p = 0:pars$pmax)
    b_vals <- slice.index(arr, 2) - 1
    y_vals <- slice.index(arr, 1) - 1
    w1 <- function(a, b) sum(abs(cumsum(a) - cumsum(b)))
    y_red <- rowSums(jr)
    w1_y <- w1(y_red, apply(arr, 1, sum))
    w1_b <- w1(rev(y_red), apply(arr, 2, sum))
    rbind(cbind(pr = pr, model = "reduced", pyr = pyr, pbr = pbr, p0r = 0,
                mean_empty = 0, summarize(jr), tv_yp = NA, w1_y = NA, w1_b = NA),
          cbind(pr = pr, model = "full", pyr = fe[["pyr"]], pbr = fe[["pbr"]],
                p0r = p0r, mean_empty = sum(arr * (n - y_vals - b_vals)),
                summarize(jf), tv_yp = sum(abs(jr - jf)) / 2,
                w1_y = w1_y, w1_b = w1_b))
}

t0 <- Sys.time()
res <- do.call(rbind, lapply(prs, one))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# differences, full minus reduced, for each P_R
num <- setdiff(names(res), c("pr", "model", "tv_yp", "w1_y", "w1_b"))
diffs <- do.call(rbind, lapply(prs, function(pr) {
    r <- res[res$pr == pr & res$model == "reduced", num]
    f <- res[res$pr == pr & res$model == "full", num]
    cbind(pr = pr, f - r,
          res[res$pr == pr & res$model == "full", c("tv_yp", "w1_y", "w1_b")])
}))

write.csv(res, "claude-checks/full_vs_reduced_output.csv", row.names = FALSE)
out <- c(sprintf("full_vs_reduced.R, %s, sweetsoursong %s, %.1f min",
                 format(Sys.time(), "%Y-%m-%d %H:%M"),
                 as.character(packageVersion("sweetsoursong")), elapsed),
         "", "== Values ==",
         capture.output(print(res, digits = 5, row.names = FALSE)),
         "", "== Full minus reduced ==",
         capture.output(print(diffs, digits = 3, row.names = FALSE)))
writeLines(out, "claude-checks/full_vs_reduced_output.txt")
cat(out, sep = "\n")
