int_check <- function(x, l, n, .min = NULL, .max = NULL) {
    if (!(length(x) %in% l && is.numeric(x) && x %% 1 == 0)) {
        stop("ERROR: ", n, "is not an integer of one of these length(s): ",
             paste(l, collapse = ", "))
    }
    if (!is.null(.min) && any(x < .min)) {
        stop("ERROR: ", n, " has value(s) less than ", .min)
    }
    if (!is.null(.max) && any(x > .max)) {
        stop("ERROR: ", n, " has value(s) greater than ", .max)
    }
    invisible(NULL)
}
dbl_check <- function(x, l, n, .min = NULL, .max = NULL) {
    if (!(length(x) %in% l && is.numeric(x))) {
        stop("ERROR: ", n, " is not a numeric of one of these length(s): ",
             paste(l, collapse = ", "))
    }
    if (!is.null(.min) && any(x < .min)) {
        stop("ERROR: ", n, " has value(s) less than ", .min)
    }
    if (!is.null(.max) && any(x > .max)) {
        stop("ERROR: ", n, " has value(s) greater than ", .max)
    }
    invisible(NULL)
}




#' Simulate plant metacommunities with stochasticity.
#'
#' @inheritParams plant_metacomm
#' @param n_sigma Sample size of the approximated binomial distribution that
#'     generates stochasticity at each time step.
#' @param season_len Length of seasons. Set to `NULL` to not have seasons.
#' @param season_surv Survival of microbes from the end of one season to
#'     the start of the next.
#' @param season_sigma Standard deviation of the normal distribution used to
#'     generate between-season variation in logit-transformed relative
#'     species compositions.
#' @param n_reps Number of reps to simulate.
#'
#' @return A data frame of yeast, bacteria, and pollinator densities at each
#'     plant through time and for each repetition.
#'
#' @export
#'
#'
plant_metacomm_stoch <- function(np,
                                 u,
                                 Y0 = NULL,
                                 B0 = NULL,
                                 m = 0.1,
                                 d_yp = 1.1,
                                 d_b0 = 0.3,
                                 d_bp = 0.4,
                                 g_yp = 0.005,
                                 g_b0 = 0.02,
                                 g_bp = 0.002,
                                 L_0 = 1 / np,
                                 X = 0,
                                 n_sigma = 200,
                                 season_len = 150,
                                 season_surv = 0.1,
                                 season_sigma = 0,
                                 n_reps = 100,
                                 dt = 0.1,
                                 max_t = 3000,
                                 closed = TRUE) {

    NotQuiteZero <- 1e-6  # used below for checking

    int_check(np, 1, "np", .min = 2)

    if ((!is.null(Y0) && is.null(B0)) || (is.null(Y0) && !is.null(B0))) {
        stop("ERROR: If Y0 is NULL, B0 must also be (and vice versa).")
    }
    if (is.null(Y0)) {
        Y0 <- seq(0.1, 0.4, length.out = np)
        B0 <- 0.5 - Y0
    }

    dbl_check(u, 1, "u", .min = 0)
    dbl_check(Y0, c(1, np), "Y0", .min = 0, .max = 1)
    dbl_check(B0, c(1, np), "B0", .min = 0, .max = 1)

    dbl_check(m, c(1, np), "m", .min = 0)
    dbl_check(d_yp, c(1, np), "d_yp", .min = 0)
    dbl_check(d_b0, c(1, np), "d_b0", .min = 0)
    dbl_check(d_bp, c(1, np), "d_bp", .min = 0)
    dbl_check(g_yp, c(1, np), "g_yp", .min = 0)
    dbl_check(g_b0, c(1, np), "g_b0", .min = 0)
    dbl_check(g_bp, c(1, np), "g_bp", .min = 0)
    dbl_check(L_0, c(1, np), "L_0", .min = 0)
    dbl_check(X, 1, "X", .min = 0)

    dbl_check(n_sigma, 1, "n_sigma", .min = 1)
    if (!is.null(season_len)) dbl_check(season_len, 1, "season_len", .min = 1)
    dbl_check(season_surv, 1, "season_surv", .min = NotQuiteZero)
    dbl_check(season_sigma, 1, "season_sigma", .min = 0)
    int_check(n_reps, 1, "n_reps", .min = 1)

    dbl_check(dt, 1, "dt", .min = NotQuiteZero)
    dbl_check(max_t, 1, "max_t", .min = 1)
    stopifnot(length(closed) == 1 && is.logical(closed))

    if (is.null(season_len)) season_len <- max_t + 1.0

    if (closed) {
        g_yp <- 0
        g_b0 <- 0
        g_bp <- 0
    }

    if (length(Y0) == 1) Y0 <- rep(Y0, np)
    if (length(B0) == 1) B0 <- rep(B0, np)
    if (length(m) == 1) m <- rep(m, np)
    if (length(d_yp) == 1) d_yp <- rep(d_yp, np)
    if (length(d_b0) == 1) d_b0 <- rep(d_b0, np)
    if (length(d_bp) == 1) d_bp <- rep(d_bp, np)
    if (length(g_yp) == 1) g_yp <- rep(g_yp, np)
    if (length(g_b0) == 1) g_b0 <- rep(g_b0, np)
    if (length(g_bp) == 1) g_bp <- rep(g_bp, np)
    if (length(L_0) == 1) L_0 <- rep(L_0, np)

    if (all(u == 0)) {
        .P <- 1 / np
        .PLP <- (.P / (L_0 + .P))
        if (all((d_yp * .PLP) == (d_b0 + d_bp * .PLP))) {
            warning("You may be at an unstable equilibrium")
        }
    }

    plant_metacomm_stoch_cpp(n_reps, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                             L_0, u, X, Y0, B0,
                             n_sigma, season_len, season_surv, season_sigma,
                             dt, max_t) |>
        tibble::as_tibble() |>
        dplyr::mutate(p = factor(as.integer(p), levels = 0:(np-1L),
                                 labels = paste("patch", 1:np)),
                      rep = factor(as.integer(rep), levels = 1:n_reps))

}
