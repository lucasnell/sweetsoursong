# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

plant_metacomm_stoch_cpp <- function(n_reps, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X, Y0, B0, n_sigma, season_len, season_surv, q, open_sys, dt, max_t, burnin, save_every, begin_end, summarize) {
    .Call(`_sweetsoursong_plant_metacomm_stoch_cpp`, n_reps, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X, Y0, B0, n_sigma, season_len, season_surv, q, open_sys, dt, max_t, burnin, save_every, begin_end, summarize)
}

plant_metacomm_cpp <- function(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X, Y0, B0, open_sys, dt, max_t) {
    .Call(`_sweetsoursong_plant_metacomm_cpp`, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X, Y0, B0, open_sys, dt, max_t)
}

