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
#' @param q Single numeric dictating to what extent new-season
#'     abundances are affected by the previous season's abundances versus
#'     random noise.
#'     Defaults to `1`.
#' @param n_reps Number of reps to simulate.
#' @param burnin Number of days to count as "burn-in" and not record
#'     in output. This can help to avoid vector memory limit error in
#'     simulations that use many plants.
#' @param save_every Output will be stored every `save_every` days.
#'     This argument is ignored if `summarize = "rep"`.
#'     Must be a multiple of `dt`.
#'     Defaults to `NULL`, which results in all time points being saved.
#' @param summarize Single string for how and whether to summarize output.
#'     If `"none"`, no summarizing happens, and the output is a time series of
#'     microbial abundances and pollinator visits by plant.
#'     If `"time"`, summarizing happens by rep and time point, so you'll get a
#'     single row of output per rep and time step.
#'     Each row will have the mean dissimilarity between plant communities
#'     (`"BC"`), mean species diversity (`"H"`), total yeast abundance
#'     across all plants (`"sumY"`), and total bacteria abundance (`"sumB"`).
#'     If `"rep"`, summarizing happens by rep, so you'll get a single row
#'     of output per rep.
#'     Each row will have the mean dissimilarity between plant communities
#'     with median taken through time (`"BC"`),
#'     mean species diversity (`"H"`), minimum total yeast abundance
#'     through time (`"minY"`), and minimum total bacteria abundance (`"minB"`).
#'     Defaults to `"none"`.
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
                                 L_0 = 0.5,
                                 X = 0,
                                 n_sigma = 200,
                                 season_len = 150,
                                 season_surv = 0.1,
                                 q = 1,
                                 n_reps = 100,
                                 dt = 0.1,
                                 max_t = 3000,
                                 burnin = 0,
                                 save_every = NULL,
                                 no_immig = TRUE,
                                 open_sys = TRUE,
                                 summarize = "none") {

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
    dbl_check(q, 1, "q", .min = 0, .max = 1)
    int_check(n_reps, 1, "n_reps", .min = 1)

    dbl_check(dt, 1, "dt", .min = NotQuiteZero)
    dbl_check(max_t, 1, "max_t", .min = 1)
    dbl_check(burnin, 1, "burnin")
    if (is.null(save_every)) save_every <- dt
    dbl_check(save_every, 1, "save_every")
    if (save_every %% dt != 0) stop("`save_every` must be a multiple of `dt`")
    stopifnot(length(no_immig) == 1 && is.logical(no_immig))
    stopifnot(length(summarize) == 1 && is.character(summarize))
    summarize <- match.arg(summarize, c("none", "time", "rep"))
    summarize <- which(c("none", "time", "rep") == summarize) - 1L

    if (summarize == 2L) save_every <- dt

    if (is.null(season_len)) season_len <- max_t + 1.0

    if (no_immig) {
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

    out_df <- plant_metacomm_stoch_cpp(n_reps, m, d_yp, d_b0, d_bp,
                                       g_yp, g_b0, g_bp,
                                       L_0, u, X, Y0, B0, n_sigma,
                                       season_len, season_surv, q,
                                       open_sys, dt, max_t, burnin, save_every,
                                       summarize) |>
        tibble::as_tibble(.name_repair = \(x) make.names(x, TRUE)) |>
        # this col will eventually be "rep"
        dplyr::mutate(X = factor(as.integer(X), levels = 1:n_reps))
    if (summarize == 0L) {  # "none"
        colnames(out_df) <- c("rep", "t", "p", "Y", "B", "P")
        out_df$p <- factor(as.integer(out_df$p), levels = 0:(np-1L),
                           labels = paste("patch", 1:np))
    } else if (summarize == 1L) {  # "time"
        colnames(out_df) <- c("rep", "t", "BC", "H", "sumY", "sumB")
    } else {  # "rep"
        colnames(out_df) <- c("rep", "BC", "H",
                              "minY", "maxY", "meanY", "meanlogY",
                              "minB", "maxB", "meanB", "meanlogB")
    }
    return(out_df)

}
