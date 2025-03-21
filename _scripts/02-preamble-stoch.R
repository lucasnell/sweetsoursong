
#'
#' Preamble for scripts used for stochastic simulations.
#'

source("_scripts/00-preamble.R")

suppressPackageStartupMessages({
    library(ggtext)
    # library(grid)
    library(RcppParallel)
})
setThreadOptions(numThreads = max(defaultNumThreads() - 2L, 1L))



# data file created in stochastic simulations:
stoch_rds_files <- list(
    # sims of outcomes / abundances:
    main = "_data/big-stoch-sims.rds",
    # sims of mutual invasibility:
    mutual_inv = "_data/big-stoch-mutual-inv-sims.rds",
    # sims of mutual invasiblity - growth rates:
    mutual_inv_growth = "_data/big-stoch-mutual-inv-growth-sims.rds")


#' "coexistence" is the boundary where coexistence occurs at u=0.
#' The other two are on either side of that boundary and produce competitive
#' exclusion.
d_yp__ = c("bacteria only" = 0.8,
           "coexistence" = 1,
           "yeast only" = 1.2)



#'
#' Convert columns into factors, optionally excluding some.
#' It also makes the species column prettier for plots.
#'
add_factors <- function(d, .exclude = NULL) {
    if ("d_yp" %in% colnames(d) && !is.factor(d[["d_yp"]]) &&
        ! "d_yp" %in% .exclude) {
        d <- mutate(d, d_yp = factor(d_yp, levels = unname(d_yp__),
                                     labels = names(d_yp__)))
    }
    if ("q" %in% colnames(d) && !is.factor(d[["q"]]) &&
        ! "q" %in% .exclude) {
        d <- mutate(d, q = factor(q))
    }
    if ("u" %in% colnames(d) && !is.factor(d[["u"]]) &&
        ! "u" %in% .exclude) {
        d <- mutate(d, u = factor(u))
    }
    adj_sp <- function(sp_col, col_name) {
        if (all(sp_col %in% c("Y", "B"))) {
            new_sp <- factor(sp_col, levels = c("Y", "B"),
                             labels = c("yeast", "bacteria"))
        } else if (all(sp_col %in% c("yeast", "bacteria"))) {
            new_sp <- factor(sp_col, levels = c("yeast", "bacteria"))
        } else if (all(sp_col %in% c("Y", "B", "P"))) {
            new_sp <- factor(sp_col, levels = c("Y", "B", "P"),
                             labels = c("yeast", "bacteria", "pollinators"))
        } else if (all(sp_col %in% c("yeast", "bacteria", "pollinators"))) {
            new_sp <- factor(sp_col, levels = c("yeast", "bacteria",
                                                "pollinators"))
        } else {
            stop("strange values in d$", col_name)
        }
        return(new_sp)
    }
    for (spc in c("species", "rare_sp", "inv_sp", "res_sp")) {
        if (spc %in% colnames(d) && !is.factor(d[[spc]]) &&
            ! spc %in% .exclude) {
            d[[spc]] <- adj_sp(d[[spc]], spc)
        }
    }
    return(d)
}




