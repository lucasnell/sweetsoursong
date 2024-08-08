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



#' Simulate plant metacommunities.
#'
#' @param np Number of plants.
#' @param u Affects how strongly pollinators are attracted to plants with the
#'     greatest proportion of palatable nectar.
#' @param Y0 Starting yeast proportions.
#' @param B0 Starting bacteria proportions.
#' @param m Flower senescence.
#' @param d_yp Pollinator-dependent, within-plant dispersal rate for yeast.
#' @param d_b0 Pollinator-independent, within-plant dispersal rate for bacteria.
#' @param d_bp Pollinator-dependent, within-plant dispersal rate for bacteria.
#' @param g_yp Pollinator-dependent rate of dispersal from sources other
#'     than the focal plant for yeast.
#'     This is overridden by the `closed` parameter by default!
#' @param g_b0 Pollinator-independent rate of dispersal from sources other
#'     than the focal plant for bacteria.
#'     This is overridden by the `closed` parameter by default!
#' @param g_bp Pollinator-dependent rate of dispersal from sources other
#'     than the focal plant for bacteria.
#'     This is overridden by the `closed` parameter by default!
#' @param L_0 Half saturation ratio for the effect of pollinators on dispersal.
#' @param X Attraction of pollinators to non-focal plants.
#' @param dt Time step for ODE solver.
#' @param max_t Number of days to simulate.
#' @param closed Single logical for a closed system, which will override
#'     all the `g_*` parameters and set them to zero. Defaults to TRUE.
#'
#' @return A data frame of yeast, bacteria, and pollinator densities at each
#'     plant through time.
#'
#' @export
#'
#'
plant_metacomm <- function(np,
                           u,
                           Y0,
                           B0,
                           m = 0.1,
                           d_yp = 1.1,
                           d_b0 = 0.3,
                           d_bp = 0.4,
                           g_yp = 0.005,
                           g_b0 = 0.02,
                           g_bp = 0.002,
                           L_0 = 1 / np,
                           X = 0,
                           dt = 0.1,
                           max_t = 500,
                           closed = TRUE) {

    int_check(np, 1, "np", .min = 2)
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
    dbl_check(dt, 1, "dt", .min = 0)
    dbl_check(max_t, 1, "max_t", .min = 0)
    stopifnot(length(closed) == 1 && is.logical(closed))

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

    plant_metacomm_cpp(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X,
                       Y0, B0, dt, max_t) |>
        tibble::as_tibble() |>
        dplyr::mutate(p = factor(p, levels = 0:(np-1L),
                                 labels = paste("patch", 1:np)))

}
