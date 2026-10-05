test_that("flower_rates match the Mathematica rates", {
    for (mode in modes) {
        ref <- load_ref(mode)
        if (is.null(ref)) next
        p <- ss_params()
        for (r in ref$rates) {
            got <- flower_rates(r$Y, r$B, r$P, r$P0R, r$PYR, r$PBR, p)
            expect_equal(as.numeric(got), num(r$rates), tolerance = 1e-14)
        }
    }
})

test_that("reduced rates sum to the reduced ODE", {
    p <- ss_params()
    r <- reduced_rates(10, 2, 30, 20, p)
    rhs <- determ_rhs(10, 2, 30, 20, p, immigration = TRUE)
    expect_equal(rhs[, "dy"],
                 unname(r[, "y_plus"] + r[, "py_plus"] - r[, "b_plus"] - r[, "pb_plus"]))
})

test_that("p_crit is 1 at the default parameters", {
    expect_equal(p_crit(), 1)
})
