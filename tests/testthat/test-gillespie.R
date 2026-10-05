test_that("Gillespie jumps follow the closure rates", {
    skip_on_cran()
    p <- ss_params()
    pyr <- 18.219389295557
    pbr <- p$n - pyr
    hz <- vacancy_hazards(pyr, pbr, p)
    set.seed(12)
    sim <- gillespie_reduced(0, 0, 2e4, pyr, pbr, p, hazards = hz)
    expect_equal(tail(sim$t, 1), 2e4)
    expect_true(all(sim$y >= 0 & sim$y <= p$n & sim$p >= 0))
    # at each y, the share of y jumps that go up should be up / (up + down)
    k <- seq_len(nrow(sim) - 2)
    dy <- diff(sim$y)[k]
    from <- sim$y[k]
    jumps <- dy != 0
    up_share <- tapply(dy[jumps] > 0, from[jumps], mean)
    n_obs <- tapply(dy[jumps] > 0, from[jumps], length)
    ys <- as.integer(names(n_obs))[n_obs > 200]
    expect_gt(length(ys), 3)
    for (y in ys) {
        obs <- up_share[[as.character(y)]]
        expd <- hz$up[y + 1] / (hz$up[y + 1] + hz$down[y + 1])
        se <- sqrt(expd * (1 - expd) / n_obs[[as.character(y)]])
        expect_lte(abs(obs - expd), 5 * se)
    }
})

test_that("Gillespie time averages match the stationary distribution", {
    skip_on_cran()
    p <- ss_params()
    hz <- vacancy_hazards(40, 10, p)
    pi_y <- birth_death_stationary(hz$up, hz$down)
    set.seed(1)
    t_end <- 4e4
    sim <- gillespie_reduced(1, 1, t_end, 40, 10, p, hazards = hz)
    # batch means over 20 equal time windows after a burn-in
    edges <- seq(2e3, t_end, length.out = 21)
    bm <- vapply(1:20, function(i) {
        s <- sim[sim$t >= edges[i] & sim$t <= edges[i + 1], ]
        temporal_mean(s)[["y"]]
    }, 0)
    se <- stats::sd(bm) / sqrt(20)
    expect_lt(abs(mean(bm) - sum(pi_y * 0:p$n)), 4 * se)
})
