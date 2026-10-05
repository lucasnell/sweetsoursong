test_that("closed-metacommunity examples (Figs 3-4) match Mathematica", {
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        p <- ss_params()
        with_mode(mode, {
            for (ex in ref$examples) {
                s <- p$n * ex$PR
                x <- meta_equilibrium(ex$PR, ex$start, p)
                # FindRoot stops within ~1e-8 of a root on a bound
                expect_lt(abs(x - ex$pyr) / s, 1e-8)
                x <- ex$pyr
                j <- joint_distribution(x, s - x, p)
                expect_lt(max(abs(j - matrix(num(ex$joint), nrow(j), byrow = TRUE))), 1e-12)
                expect_equal(unname(joint_means(j)), num(ex$mean), tolerance = 1e-10)
                d <- decompose_distribution(j)
                # Without Chop, single-cell minima in underflowing tails can
                # differ by round-off (PR = 0.6); labels still agree.
                if (mode == "chop") {
                    mk <- regional_minima(-j)
                    expect_equal(as.vector(mk > 0),
                                 as.vector(matrix(num(ex$markers), nrow(mk), byrow = TRUE) > 0))
                }
                expect_equal(as.vector(d$labels),
                             as.vector(matrix(num(ex$labels), nrow(d$labels), byrow = TRUE)))
                expect_lt(max(abs(d$weights - num(ex$weights))), 1e-10)
                cm <- matrix(num(ex$compMeans), ncol = 2, byrow = TRUE)
                expect_lt(max(abs(d$means - cm)), 1e-9)
            }
        })
    }
})

test_that("Fig 5 continuation matches Mathematica", {
    skip_on_cran()
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        p <- ss_params()
        f5 <- ref$fig5
        prs <- vapply(f5, `[[`, 0, "PR")
        pyr <- vapply(f5, `[[`, 0, "pyr")
        with_mode(mode, {
            x1 <- meta_equilibrium(1, 25, p)
            run <- function(grid, x) vapply(grid, function(pr) {
                x <<- meta_equilibrium(pr, x, p)
                x
            }, 0)
            dn <- seq(0.99, 0.74, by = -0.01)
            up <- seq(1, 2.89, by = 0.01)
            got <- c(run(dn, x1), run(up, x1))
            got <- got[order(c(dn, up))]
            expect_lt(max(abs(got - pyr) / pmax(1, pyr)), 1e-9)
            for (k in c(1, 50, 100, 150, 216)) {
                j <- joint_distribution(pyr[k], p$n * prs[k] - pyr[k], p)
                d <- decompose_distribution(j)
                expect_lt(max(abs(d$weights - num(f5[[k]]$weights))), 1e-9)
            }
        })
    }
})
