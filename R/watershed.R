# Neighbour index matrix for a d-dimensional grid: one row per cell, one
# column per neighbour (8 in 2-D, 26 in 3-D, diagonals included), NA outside
# the grid. Offsets are ordered with the first dimension varying fastest.
neighbour_matrix <- function(dims) {
    d <- length(dims)
    off <- as.matrix(expand.grid(rep(list(-1:1), d)))
    off <- off[rowSums(abs(off)) > 0, , drop = FALSE]
    sub <- arrayInd(seq_len(prod(dims)), dims)
    mult <- c(1, cumprod(dims)[-d])
    nb <- matrix(NA_integer_, nrow(sub), nrow(off))
    for (k in seq_len(nrow(off))) {
        s <- sweep(sub, 2, off[k, ], `+`)
        ok <- rowSums(s < 1 | sweep(s, 2, dims, `>`)) == 0
        nb[ok, k] <- as.integer(1 + (s[ok, , drop = FALSE] - 1) %*% mult)
    }
    nb
}

# Cells in row-major ("raster") order: last index varies fastest, as an
# image is scanned row by row.
raster_order <- function(dims) {
    d <- length(dims)
    grid <- as.matrix(expand.grid(rev(lapply(dims, seq_len))))[, d:1, drop = FALSE]
    as.integer(1 + (grid - 1) %*% c(1, cumprod(dims)[-d]))
}

# Union-find root with path halving.
uf_find <- function(parent, i) {
    while (parent[i] != i) i <- parent[i]
    i
}

#' Regional minima of an array
#'
#' A regional minimum is a connected set of cells with equal value whose
#' neighbours outside the set all have larger values (Wolfram `MinDetect`
#' with h = 0). Neighbours include diagonals.
#'
#' @param x Numeric matrix or 3-D array.
#' @return Integer array of the same shape: 0 outside minima, and minima
#'   labelled 1, 2, ... in raster order of their first cell.
#' @export
regional_minima <- function(x) {
    dims <- dim(x)
    nb <- neighbour_matrix(dims)
    xv <- as.numeric(x)
    nv <- matrix(xv[nb], nrow(nb))
    has_lower <- rowSums(nv < xv, na.rm = TRUE) > 0
    # group equal-valued neighbouring cells into plateaus
    parent <- seq_along(xv)
    eq <- which(!is.na(nv) & nv == xv, arr.ind = TRUE)
    if (nrow(eq) > 0) {
        a <- eq[, 1]
        b <- nb[eq]
        for (k in seq_along(a)) {
            ra <- uf_find(parent, a[k])
            rb <- uf_find(parent, b[k])
            if (ra != rb) parent[max(ra, rb)] <- min(ra, rb)
        }
    }
    root <- vapply(seq_along(xv), function(i) uf_find(parent, i), 1L)
    plateau_low <- tapply(has_lower, root, any)
    is_min <- !plateau_low[as.character(root)]
    lab <- integer(length(xv))
    k <- 0L
    seen <- integer(0)
    for (i in raster_order(dims)) {
        if (!is_min[i]) next
        r <- root[i]
        if (!(r %in% seen)) {
            k <- k + 1L
            seen <- c(seen, r)
            lab[root == r] <- k
        }
    }
    array(lab, dims)
}

#' Watershed basins
#'
#' Each cell follows its steepest descent (the neighbour with the most
#' negative difference in value, diagonals included) until it reaches a
#' labelled minimum, and takes that minimum's label. This reproduces Wolfram
#' `WatershedComponents[image, markers, Method -> "Basins"]` on the
#' distributions in this package. Cells on non-minimal plateaus, which have
#' no lower neighbour, take the label of a neighbour.
#'
#' @param x Numeric matrix or 3-D array (minima are basin bottoms).
#' @param markers Integer array from [regional_minima()].
#' @return Integer array of basin labels.
#' @export
watershed_basins <- function(x, markers = regional_minima(x)) {
    dims <- dim(x)
    nb <- neighbour_matrix(dims)
    xv <- as.numeric(x)
    nv <- matrix(xv[nb], nrow(nb))
    g <- nv - xv
    g[is.na(g)] <- Inf
    k <- max.col(-g, ties.method = "first")
    step_to <- nb[cbind(seq_along(xv), k)]
    downhill <- g[cbind(seq_along(xv), k)] < 0
    nxt <- ifelse(downhill, step_to, seq_along(xv))
    lab <- as.integer(markers)
    # follow descent paths; repeat until labels stop changing
    repeat {
        new <- ifelse(lab > 0L, lab, lab[nxt])
        if (identical(new, lab)) break
        lab <- new
    }
    # cells that end on a non-minimal plateau: take a labelled neighbour's
    # label, spreading until all are labelled
    while (any(lab == 0L)) {
        z <- which(lab == 0L)
        nl <- matrix(lab[nb[z, , drop = FALSE]], length(z))
        nl[is.na(nl)] <- 0L
        got <- apply(nl, 1, function(v) if (any(v > 0)) v[v > 0][1] else 0L)
        if (all(got == 0L)) break
        lab[z] <- got
    }
    array(lab, dims)
}

#' Split a bimodal distribution into its modes
#'
#' Assigns each cell of a joint distribution to the basin of one mode by a
#' watershed on the negative probabilities (Chris Klausmeier's
#' `DecomposeDistribution`), and returns each mode's probability mass and
#' within-mode means.
#'
#' @param joint Matrix from [joint_distribution()] (rows y, columns p) or a
#'   3-D array from [full_distribution()] (y, b, p). Dimnames give the state
#'   values.
#' @param floor Probabilities below `floor` are treated as 0 when locating
#'   modes (weights and means still use the full distribution). Without
#'   Chop, round-off of order 1e-17 in the near-empty tails of the full
#'   model creates spurious minima; `floor = 1e-10` removes them, matching
#'   the effect of Chop in Chris's code. Default 0 (no floor).
#' @return A list with `weights` (mass of each mode), `means` (one row per
#'   mode), and `labels` (basin of each cell). Modes are in the order of
#'   Chris's code: decreasing label, so the mode found last in raster order
#'   comes first.
#' @export
decompose_distribution <- function(joint, floor = 0) {
    dat <- joint
    ws <- dat
    ws[ws < floor] <- 0
    labels <- watershed_basins(-ws)
    ncomp <- max(labels)
    comps <- rev(seq_len(ncomp))
    tot <- sum(dat)
    vals <- lapply(seq_along(dim(dat)), function(k) {
        v <- as.numeric(dimnames(dat)[[k]])
        array(v[slice.index(dat, k)], dim(dat))
    })
    weights <- vapply(comps, function(k) sum(dat[labels == k]) / tot, 0)
    means <- t(vapply(comps, function(k) {
        w <- dat[labels == k]
        vapply(vals, function(g) sum(g[labels == k] * w) / sum(w), 0)
    }, numeric(length(vals))))
    colnames(means) <- names(dimnames(dat))
    list(weights = weights, means = means, labels = labels)
}
