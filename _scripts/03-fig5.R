# Fig 5 (panels B-F): community structure along P_R in the closed
# metacommunity. Panel A is a schematic. Continuation as in section 2.3 of
# exploration_unified_10x.nb: from the root at P_R = 1, up to 2.89 and down
# to 0.74 in steps of 0.01, each step starting from the previous root.

source("_scripts/00-preamble.R")

res <- cached("fig5-continuation.rds", {
    x1 <- meta_equilibrium(1, 25, pars)
    one <- function(pr, x) {
        pyr <- meta_equilibrium(pr, x, pars)
        pbr <- pars$n * pr - pyr
        j <- joint_distribution(pyr, pbr, pars)
        d <- decompose_distribution(j)
        # the first mode is the yeast-dominated one when there are two
        iy <- which.max(d$means[, "y"])
        ib <- which.min(d$means[, "y"])
        data.frame(pr = pr, pyr = pyr,
                   mean_y = joint_means(j)[["y"]], mean_p = joint_means(j)[["p"]],
                   n_modes = length(d$weights),
                   p_y = d$weights[iy], p_b = if (ib != iy) d$weights[ib] else 0,
                   y_y = d$means[iy, "y"], p_on_y = d$means[iy, "p"],
                   y_b = d$means[ib, "y"], p_on_b = d$means[ib, "p"],
                   inv_y = inv_y(pr, pars), inv_b = inv_b(pr, pars))
    }
    run <- function(prs, x) {
        out <- vector("list", length(prs))
        for (i in seq_along(prs)) {
            out[[i]] <- one(prs[i], x)
            x <- out[[i]]$pyr
        }
        do.call(rbind, out)
    }
    up <- run(seq(1, 2.89, by = 0.01), x1)
    dn <- run(seq(0.99, 0.74, by = -0.01), x1)
    r <- rbind(dn, up)
    r[order(r$pr), ]
})
write.csv(res, data_path("fig5-continuation.csv"), row.names = FALSE)
cat("Coexistence (both modes present) for P_R in",
    paste(range(res$pr[res$n_modes == 2]), collapse = " to "), "\n")

n <- pars$n
lo <- min(res$pr); hi <- max(res$pr)
# outside the coexistence range one species is absent
ext_b <- data.frame(pr = c(0, lo))
ext_y <- data.frame(pr = c(hi, 3.6))

pB <- ggplot(res, aes(pr, mean_p)) +
    geom_line(color = col_p) +
    geom_line(data = ext_b, aes(pr, pr), color = col_p) +
    geom_line(data = ext_y, aes(pr, pr), color = col_p) +
    labs(x = NULL, y = "Pollinators", title = "B")
pC <- ggplot(res, aes(pr)) +
    geom_line(aes(y = mean_y), color = col_y) +
    geom_line(aes(y = n - mean_y), color = col_b) +
    geom_line(data = ext_b, aes(pr, 0), color = col_y) +
    geom_line(data = ext_b, aes(pr, n), color = col_b) +
    geom_line(data = ext_y, aes(pr, n), color = col_y) +
    geom_line(data = ext_y, aes(pr, 0), color = col_b) +
    labs(x = NULL, y = "Flowers", title = "C")
pD <- ggplot(res, aes(pr)) +
    geom_line(aes(y = p_y), color = col_y) +
    geom_line(aes(y = p_b), color = col_b, linetype = 2) +
    geom_line(data = ext_b, aes(pr, 0), color = col_y) +
    geom_line(data = ext_b, aes(pr, 1), color = col_b, linetype = 2) +
    geom_line(data = ext_y, aes(pr, 1), color = col_y) +
    geom_line(data = ext_y, aes(pr, 0), color = col_b, linetype = 2) +
    labs(x = NULL, y = "Proportion", title = "D")
pE <- ggplot(res, aes(pr)) +
    geom_line(aes(y = p_on_y), color = col_p) +
    geom_line(aes(y = p_on_b), color = col_p, linetype = 2) +
    geom_line(data = ext_b, aes(pr, pr), color = col_p, linetype = 2) +
    geom_line(data = ext_y, aes(pr, pr), color = col_p) +
    labs(x = NULL, y = "Pollinators", title = "E")
pF <- ggplot(res, aes(pr)) +
    geom_line(aes(y = y_y), color = col_y) +
    geom_line(aes(y = n - y_y), color = col_b) +
    geom_line(aes(y = y_b), color = col_y, linetype = 2) +
    geom_line(aes(y = n - y_b), color = col_b, linetype = 2) +
    geom_line(data = ext_b, aes(pr, 0), color = col_y, linetype = 2) +
    geom_line(data = ext_b, aes(pr, n), color = col_b, linetype = 2) +
    geom_line(data = ext_y, aes(pr, n), color = col_y) +
    geom_line(data = ext_y, aes(pr, 0), color = col_b) +
    labs(x = "Regional pollinator abundance (P_R)", y = "Flowers", title = "F")
fig5 <- (pB / pC / pD / pE / pF) & coord_cartesian(xlim = c(0, 3.6))
ggsave(fig_path("fig5-meta-PR.pdf"), fig5, width = 3.5, height = 9)
