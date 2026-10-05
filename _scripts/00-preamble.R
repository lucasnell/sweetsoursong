# Shared setup for the figure scripts. Run scripts from the project root,
# e.g. `Rscript _scripts/01-fig2.R`. Figures go to `_figures/` and cached
# results to `_data/` (both git-ignored).

suppressPackageStartupMessages({
    library(sweetsoursong)
    library(ggplot2)
    library(patchwork)
})

# No Chop (package default). Figs 2-6 are identical with Chop on
# (claude-checks/r-port/compare_chop_output.txt); set TRUE to reproduce
# Chris Klausmeier's Mathematica code exactly.
options(sweetsoursong.chop = FALSE)

pars <- ss_params()

n_cores <- max(1L, parallel::detectCores() - 2L)

col_y <- rgb(0.880722, 0.611041, 0.142051)
col_b <- rgb(0.368417, 0.506779, 0.709798)
col_p <- rgb(0.5, 0, 0.5)

theme_set(theme_classic(base_size = 10) +
              theme(strip.background = element_blank(),
                    plot.title = element_text(size = 10)))

fig_path <- function(name) file.path("_figures", name)
data_path <- function(name) file.path("_data", name)

# Cache an expensive result as .rds.
cached <- function(name, expr) {
    path <- data_path(name)
    if (file.exists(path)) return(readRDS(path))
    out <- force(expr)
    saveRDS(out, path)
    out
}

# Nullclines of the deterministic one-plant model in the (P, Y) plane.
# Y nullcline: contour of dY/dt = 0 on a grid; P nullcline: P = E[P | Y].
nullclines <- function(pyr, pbr, pars, immigration, p_max = 3,
                       ny = 401, np = 401) {
    ys <- seq(0, pars$n, length.out = ny)
    ps <- seq(1e-6, p_max, length.out = np)
    g <- expand.grid(p = ps, y = ys)
    dy <- determ_rhs(g$y, g$p, pyr, pbr, pars, immigration)[, "dy"]
    z <- matrix(dy, np, ny)
    cl <- grDevices::contourLines(ps, ys, z, levels = 0)
    ync <- do.call(rbind, lapply(seq_along(cl), function(i) {
        data.frame(p = cl[[i]]$x, y = cl[[i]]$y, piece = i)
    }))
    pnc <- data.frame(y = ys, p = p_mean_given_y(ys, pyr, pbr, pars))
    pnc <- pnc[pnc$p <= p_max, ]
    list(y = ync, p = pnc)
}

# Phase-plane panel: nullclines and equilibria, optionally over a joint
# stationary distribution (grayscale).
phase_panel <- function(pyr, pbr, pars, immigration, p_max = 3,
                        joint = NULL) {
    nc <- nullclines(pyr, pbr, pars, immigration, p_max)
    if (!immigration) {
        # without immigration, Y = 0 and Y = n are Y-nullclines everywhere
        edge <- data.frame(p = rep(c(0, p_max), 2), y = rep(c(0, pars$n), each = 2),
                           piece = rep(c(-1, -2), each = 2))
        nc$y <- rbind(nc$y, edge)
    }
    eq <- determ_equilibria(pyr, pbr, pars, immigration)
    n <- pars$n
    gg <- ggplot()
    if (!is.null(joint)) {
        jd <- expand.grid(y = 0:pars$n, p = 0:pars$pmax)
        jd$prob <- as.numeric(joint)
        jd <- jd[jd$p <= p_max + 0.5, ]
        gg <- gg + geom_tile(data = jd, aes(p, y, fill = prob)) +
            scale_fill_gradient(low = "white", high = "black", guide = "none")
    }
    gg +
        geom_path(data = nc$y, aes(p, y, group = piece), color = col_y) +
        geom_path(data = nc$p, aes(p, y), color = col_p) +
        geom_point(data = eq, aes(p, y, shape = stable), size = 1.8,
                   fill = "white") +
        scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
        scale_y_continuous("Yeast (Y)",
                           breaks = seq(0, n, 10),
                           sec.axis = sec_axis(transform = identity,
                                               breaks = seq(0, n, 10),
                                               labels = function(b) n - b,
                                               name = "Bacteria (B)")) +
        scale_x_continuous("Pollinators (P)") +
        coord_cartesian(xlim = c(0, p_max), ylim = c(-1, n + 1), expand = FALSE) +
        theme(axis.title.y.left = element_text(color = col_y),
              axis.title.y.right = element_text(color = col_b),
              axis.title.x = element_text(color = col_p))
}

# The closed-metacommunity examples of Figs 3-4: P_R and FindRoot starting
# values from section 2.2 of Chris's notebook.
examples <- data.frame(pr = c(0.6, 1, 1.8, 2.6, 3),
                       start = c(25, 25, 85, 130, 130),
                       label = c("bacteria win", "coexistence", "coexistence",
                                 "coexistence", "yeast win"))
