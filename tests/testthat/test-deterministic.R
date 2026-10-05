test_that("Fig 2 equilibria and stability match EcoEvo", {
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        p <- ss_params()
        for (f in ref$fig2) {
            pp <- update_params(p, m_b = f$mB)
            eq <- determ_equilibria(90, 10, pp, immigration = f$model == "reduced")
            expect_equal(eq$y, num(f$eqY), tolerance = 1e-9)
            expect_equal(eq$p, num(f$eqP), tolerance = 1e-9)
            expect_equal(eq$stable, as.logical(unlist(f$stable)))
        }
    }
})

test_that("determ_sim approaches a stable equilibrium", {
    p <- update_params(ss_params(), m_b = 0.05)
    out <- determ_sim(45, 2, c(0, 2000), 90, 10, p, immigration = TRUE)
    eq <- determ_equilibria(90, 10, p, immigration = TRUE)
    hi <- eq[which.max(eq$y), ]
    expect_equal(out$y[2], hi$y, tolerance = 1e-4)
})
