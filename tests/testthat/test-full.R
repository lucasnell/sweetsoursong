test_that("full-model invasion criteria match Mathematica", {
    skip_on_cran()
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        p <- ss_params()
        with_mode(mode, {
            for (v in ref$full[c(1, 2)]) {
                # FindRoot's default precision (~8 digits) on the resident
                # equilibrium limits agreement
                expect_lt(rel_diff(full_inv_y(v$PR, p), v$invY), 1e-6)
                expect_lt(rel_diff(full_inv_b(v$PR, p), v$invB), 1e-6)
            }
        })
    }
})

test_that("full-model equilibrium at PR = 1 matches Mathematica", {
    skip_on_cran()
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        f <- ref$full_example
        p <- ss_params()
        with_mode(mode, {
            x <- full_meta_equilibrium(f$PR, c(18, 31), p)
            expect_equal(unname(x), c(f$pyr, f$pbr), tolerance = 1e-6)
            arr <- full_distribution(p$n * f$PR - f$pyr - f$pbr, f$pyr, f$pbr, p)
            g <- sweetsoursong:::full_grids(arr)
            expect_equal(c(sum(arr * g$y), sum(arr * g$b), sum(arr * g$p)),
                         num(f$mean), tolerance = 1e-8)
            if (mode == "chop") {
                d <- decompose_distribution(arr)
                expect_equal(d$weights, num(f$weights), tolerance = 1e-8)
            }
        })
    }
})

test_that("tm_stationary solves TM v = 0", {
    p <- ss_params()
    tm <- full_tm(10, 20, 20, p)
    old <- options(sweetsoursong.chop = FALSE)
    on.exit(options(old))
    v <- sweetsoursong:::tm_stationary(tm)
    expect_equal(sum(v), 1)
    expect_lt(max(abs(as.numeric(tm %*% v))), 1e-12)
    expect_equal(as.numeric(Matrix::colSums(tm)), rep(0, nrow(tm)), tolerance = 1e-9)
})

test_that("full-model invasion thresholds match Mathematica", {
    skip_on_cran()
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        th <- ref$full_thresholds
        p <- ss_params()
        with_mode(mode, {
            ty <- stats::uniroot(function(pr) full_inv_y(pr, p) - 1, c(0.7, 0.76),
                                 tol = 1e-9)$root
            tb <- stats::uniroot(function(pr) full_inv_b(pr, p) - 1, c(2.8, 2.9),
                                 tol = 1e-9)$root
            expect_equal(ty, th$invY, tolerance = 1e-6)
            expect_equal(tb, th$invB, tolerance = 1e-6)
        })
    }
})
