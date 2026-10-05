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
