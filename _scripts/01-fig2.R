# Fig 2 (panels B-E): deterministic one-plant model. Panel A is a schematic.
# Settings from section 1.4 of exploration_unified_10x.nb: P_R = 2,
# P_YR : P_BR = 9 : 1.

source("_scripts/00-preamble.R")

pyr <- 2 * pars$n * 0.9
pbr <- 2 * pars$n * 0.1

panels <- data.frame(
    m_b = c(0.01, 0.015, 0.05, 0.05),
    immigration = c(FALSE, FALSE, FALSE, TRUE),
    title = c("B  No pollinator preference (m_B = m)",
              "C  Weak pollinator preference (m_B > m)",
              "D  Strong pollinator preference (m_B >> m)",
              "E  Strong preference + microbe immigration")
)

plots <- lapply(seq_len(nrow(panels)), function(i) {
    p <- update_params(pars, m_b = panels$m_b[i])
    phase_panel(pyr, pbr, p, panels$immigration[i], p_max = 3) +
        geom_vline(xintercept = p_crit(p), linetype = 3, color = "gray50") +
        ggtitle(panels$title[i])
})

eq_table <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
    p <- update_params(pars, m_b = panels$m_b[i])
    cbind(panel = LETTERS[i + 1], m_b = panels$m_b[i],
          immigration = panels$immigration[i],
          determ_equilibria(pyr, pbr, p, panels$immigration[i]))
}))
write.csv(eq_table, data_path("fig2-equilibria.csv"), row.names = FALSE)
print(eq_table[, 1:6])

fig <- wrap_plots(plots, ncol = 2)
ggsave(fig_path("fig2-1plant-determ.pdf"), fig, width = 7, height = 6.5)
