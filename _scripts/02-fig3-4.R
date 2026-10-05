# Figs 3 and 4: one plant in the closed metacommunity (reduced model,
# vacancy closure). Settings from section 2.2 of exploration_unified_10x.nb.

source("_scripts/00-preamble.R")

# ---- closed-metacommunity equilibria at the example P_R values ----
eqs <- cached("fig3-4-equilibria.rds", {
    lapply(seq_len(nrow(examples)), function(i) {
        pr <- examples$pr[i]
        pyr <- meta_equilibrium(pr, examples$start[i], pars)
        pbr <- pars$n * pr - pyr
        joint <- joint_distribution(pyr, pbr, pars)
        list(pr = pr, pyr = pyr, pbr = pbr, joint = joint,
             decomp = decompose_distribution(joint),
             inv = c(y = inv_y(pr, pars), b = inv_b(pr, pars)))
    })
})

summ <- do.call(rbind, lapply(eqs, function(e) {
    data.frame(pr = e$pr, pyr = e$pyr, pbr = e$pbr,
               mean_y = joint_means(e$joint)[["y"]],
               mean_p = joint_means(e$joint)[["p"]],
               n_modes = length(e$decomp$weights),
               weights = paste(signif(e$decomp$weights, 4), collapse = " / "),
               inv_y = e$inv[["y"]], inv_b = e$inv[["b"]])
}))
write.csv(summ, data_path("fig3-4-equilibria.csv"), row.names = FALSE)
print(summ)

# ---- Fig 3: P_R = 1.8 ----
e18 <- eqs[[which(examples$pr == 1.8)]]
hz <- vacancy_hazards(e18$pyr, e18$pbr, pars)
sims <- cached("fig3-gillespie.rds", {
    set.seed(15)
    list(A = gillespie_reduced(49, 3, 1000, e18$pyr, e18$pbr, pars, hazards = hz),
         B = gillespie_reduced(5, 1, 1000, e18$pyr, e18$pbr, pars, hazards = hz))
})

# step functions sampled on a fine grid for plotting
sample_path <- function(sim, times = seq(0, 1000, by = 1)) {
    k <- findInterval(times, sim$t)
    data.frame(t = times, y = sim$y[k], b = sim$b[k], p = sim$p[k])
}
ts_panel <- function(sim, title) {
    d <- sample_path(sim)
    flowers <- ggplot(d, aes(t)) +
        geom_line(aes(y = y), color = col_y, linewidth = 0.3) +
        geom_line(aes(y = b), color = col_b, linewidth = 0.3) +
        scale_y_continuous("Flowers", limits = c(0, pars$n)) +
        labs(x = NULL, title = title)
    polls <- ggplot(d, aes(t, p)) +
        geom_col(fill = col_p, width = 1) +
        labs(x = "Time (days)", y = "Pollinators (P)")
    flowers / polls
}
fig3 <- (ts_panel(sims$A, "A  Start Y = 49, P = 3") |
             ts_panel(sims$B, "B  Start Y = 5, P = 1") |
             (phase_panel(e18$pyr, e18$pbr, pars, immigration = TRUE,
                          p_max = 8.5, joint = e18$joint) + ggtitle("C"))) +
    plot_layout(widths = c(1, 1, 1.2))
ggsave(fig_path("fig3-1plant-stoch.pdf"), fig3, width = 9, height = 3.5)

# ---- Fig 4: 50 random plants and the stationary distribution ----
# As in the notebook: round(50 * weight) plants from each mode, each drawn
# from that mode's part of the joint distribution.
draw_plants <- function(e, n_plants = 50) {
    d <- e$decomp
    states <- expand.grid(y = 0:pars$n, p = 0:pars$pmax)
    out <- lapply(seq_along(d$weights), function(k) {
        lab <- rev(seq_along(d$weights))[k]
        w <- as.numeric(e$joint) * (as.numeric(d$labels) == lab)
        m <- round(n_plants * d$weights[k])
        if (m == 0 || sum(w) == 0) return(NULL)
        states[sample.int(nrow(states), m, replace = TRUE, prob = w), ]
    })
    out <- do.call(rbind, out)
    out$plant <- seq_len(nrow(out))
    out
}
set.seed(777)
plants <- lapply(eqs, draw_plants)

row_panel <- function(i) {
    e <- eqs[[i]]
    pl <- plants[[i]]
    title <- sprintf("%s  P_R = %.1f: %s", LETTERS[i], e$pr, examples$label[i])
    top <- ggplot(pl, aes(plant, p)) +
        geom_col(fill = col_p, width = 0.9) +
        scale_y_continuous(NULL, limits = c(0, 8), breaks = c(0, 4, 8)) +
        scale_x_continuous(NULL, breaks = NULL) +
        ggtitle(title)
    fl <- rbind(data.frame(plant = pl$plant, type = "B", n = pars$n - pl$y),
                data.frame(plant = pl$plant, type = "Y", n = pl$y))
    bottom <- ggplot(fl, aes(plant, n, fill = type)) +
        geom_col(width = 0.9) +
        scale_fill_manual(values = c(Y = col_y, B = col_b), guide = "none") +
        scale_y_continuous(NULL, limits = c(0, pars$n)) +
        scale_x_continuous(NULL, breaks = NULL)
    right <- phase_panel(e$pyr, e$pbr, pars, immigration = TRUE, p_max = 8.5,
                         joint = e$joint)
    ((top / bottom) + plot_layout(heights = c(1, 2.5)) | right) +
        plot_layout(widths = c(3, 1))
}
fig4 <- wrap_plots(lapply(seq_along(eqs), row_panel), ncol = 1)
ggsave(fig_path("fig4-meta-equil.pdf"), fig4, width = 8, height = 12)
