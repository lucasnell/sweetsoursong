test_that("regional minima and basins on a two-well surface", {
    x <- outer(1:9, 1:5, function(i, j) pmin((i - 2)^2, (i - 8)^2 + 0.5) + 0.1 * j)
    mk <- regional_minima(x)
    expect_equal(max(mk), 2L)
    lab <- watershed_basins(x, mk)
    expect_true(all(lab[1:4, ] == 1L))
    expect_true(all(lab[6:9, ] == 2L))
})

test_that("decompose_distribution weights sum to 1", {
    p <- ss_params()
    j <- joint_distribution(18.219389295557, 50 - 18.219389295557, p)
    d <- decompose_distribution(j)
    expect_equal(sum(d$weights), 1)
    expect_equal(nrow(d$means), length(d$weights))
})
