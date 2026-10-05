test_that("vacancy chain and fill probabilities match Mathematica", {
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        p <- ss_params()
        with_mode(mode, {
            for (v in ref$vacancy) {
                sys <- vacancy_system(v$y, v$b, v$pyr, v$pbr, p)
                expect_equal(sys$y_repl, num(v$yRepl), tolerance = 1e-13)
                expect_equal(sys$b_repl, num(v$bRepl), tolerance = 1e-13)
                expect_equal(sys$p_down, num(v$pDown), tolerance = 1e-13)
                expect_equal(sys$p_up, num(v$pUp), tolerance = 1e-13)
                h <- fill_probs(v$y, v$b, v$pyr, v$pbr, p)
                expect_lt(max(abs(h[, "h_y"] - num(v$hY))), 1e-12)
                expect_lt(max(abs(h[, "h_b"] - num(v$hB))), 1e-12)
            }
        })
    }
})

test_that("hazards, pi(y), and E[YP] match Mathematica", {
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        p <- ss_params()
        with_mode(mode, {
            for (h in ref$hazards) {
                pbr <- p$n * h$PR - h$pyr
                hz <- vacancy_hazards(h$pyr, pbr, p)
                expect_lt(rel_diff(hz$up, num(h$up), 1e-300), 1e-10)
                expect_lt(rel_diff(hz$down, num(h$down), 1e-300), 1e-10)
                expect_lt(max(abs(hz$tail_mass - num(h$tail))), 1e-15)
                pi_y <- y_distribution(h$pyr, pbr, p)
                expect_lt(max(abs(pi_y - num(h$piY))), 1e-12)
                expect_lt(rel_diff(e_yp(h$pyr, h$PR, p), h$eYP, 1e-300), 1e-10)
            }
        })
    }
})

test_that("an exactly zero pool gives the one-species limit", {
    p <- ss_params()
    pr <- 3
    # yeast only: pbr = 0 exactly puts all mass at y = n
    pi_y <- y_distribution(p$n * pr, 0, p)
    expect_equal(pi_y[p$n + 1], 1)
    # close to the limit, as Mathematica's FindRoot result
    pi_near <- y_distribution(149.9999991669740, p$n * pr - 149.9999991669740, p)
    expect_gt(pi_near[p$n + 1], 0.99)
    # bacteria only: pyr = 0 exactly puts all mass at y = 0
    expect_equal(y_distribution(0, p$n * 0.6, p)[1], 1)
})
