# Fig 6: outcome of competition over P_R and m_B, from invasion criteria.
# Ranges from section 2.4 of exploration_unified_10x.nb: P_R from 0.1 to 10
# (log scale); m_B from m to 0.1 (A) and from m to 1 (B).

source("_scripts/00-preamble.R")

prs <- 10^seq(-1, 1, length.out = 121)
grid_a <- cached("fig6-grid-a.rds",
                 outcome_grid(prs, seq(pars$m, 0.1, length.out = 46), pars,
                              cores = n_cores))
grid_b <- cached("fig6-grid-b.rds",
                 outcome_grid(prs, seq(pars$m, 1, length.out = 100), pars,
                              cores = n_cores))
for (g in list(grid_a, grid_b)) print(table(g$outcome))

pc <- p_crit(pars)
ase <- data.frame(pr = prs)
ase$m_b <- pars$m * pmax(ase$pr / pc, pc / ase$pr)

outcome_cols <- c(bacteria = col_b,
                  coexist = grDevices::adjustcolor(
                      grDevices::colorRampPalette(c(col_y, col_b))(3)[2], 0.5),
                  yeast = col_y, neither = "white")

panel <- function(g, ymax, title) {
    ggplot(g, aes(pr, m_b, fill = outcome)) +
        geom_tile(height = diff(range(g$m_b)) / (length(unique(g$m_b)) - 1),
                  width = 2 / 120) +
        geom_line(data = ase[ase$m_b <= ymax, ], aes(pr, m_b), inherit.aes = FALSE,
                  linetype = 2) +
        geom_hline(yintercept = pars$m) +
        scale_fill_manual(values = outcome_cols, drop = FALSE) +
        scale_x_log10("Regional pollinator abundance (P_R)",
                      breaks = c(0.1, 0.5, 1, 5, 10)) +
        scale_y_continuous("Emigration probability, bacteria visit (m_B)") +
        coord_cartesian(ylim = c(0, ymax), expand = FALSE) +
        ggtitle(title)
}
fig6 <- panel(grid_a, 0.1, "A") / panel(grid_b, 1, "B") +
    plot_layout(guides = "collect")
ggsave(fig_path("fig6-meta-coexist.pdf"), fig6, width = 4.5, height = 7)
