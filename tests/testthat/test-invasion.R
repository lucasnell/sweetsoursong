test_that("invasion criteria match Mathematica", {
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        p <- ss_params()
        with_mode(mode, {
            for (v in ref$inv) {
                expect_lt(rel_diff(inv_y(v$PR, p), v$invY), 1e-9)
                expect_lt(rel_diff(inv_b(v$PR, p, "chris"), v$invB), 1e-9)
                expect_lt(rel_diff(inv_b(v$PR, p, "direct"), v$eYPdirectB / p$eps), 1e-9)
            }
            pm <- update_params(p, m_b = p$m)
            for (v in ref$inv_mBeqm) {
                expect_lt(rel_diff(inv_y(v$PR, pm), v$invY), 1e-9)
                # (n * PR - E[YP]) / eps cancels: round-off in E[YP] is
                # amplified by n * PR / eps
                expect_lt(abs(inv_b(v$PR, pm, "chris") - v$invB),
                          1e-12 * pm$n * v$PR / pm$eps)
            }
        })
    }
})

test_that("the two forms of inv_b are related by the m / m_b identity", {
    p <- ss_params()
    for (pr in c(0.6, 2.8)) {
        chris <- inv_b(pr, p, "chris")
        direct <- inv_b(pr, p, "direct")
        expect_equal(direct - 1, (p$m / p$m_b) * (chris - 1), tolerance = 1e-6)
    }
})

test_that("invasion thresholds match Mathematica", {
    skip_on_cran()
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        th <- ref$thresholds
        p <- ss_params()
        pm <- update_params(p, m_b = p$m)
        with_mode(mode, {
            expect_equal(inv_threshold("y", c(0.5, 1), p), th$invY_mB05, tolerance = 1e-8)
            expect_equal(inv_threshold("b", c(2.5, 3.5), p), th$invB_mB05, tolerance = 1e-8)
            expect_equal(inv_threshold("b", c(2.5, 3.5), p, "chris"), th$invB_mB05,
                         tolerance = 1e-8)
            expect_equal(inv_threshold("y", c(1.4, 1.5), pm), th$invY_mB01, tolerance = 1e-8)
            expect_equal(inv_threshold("b", c(1.5, 1.6), pm), th$invB_mB01, tolerance = 1e-8)
        })
    }
})

test_that("classify_outcome labels all four outcomes", {
    out <- classify_outcome(c(0.5, 2, 2, 0.5), c(2, 2, 0.5, 0.5))
    expect_equal(as.character(out), c("bacteria", "coexist", "yeast", "neither"))
})
