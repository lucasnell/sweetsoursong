#' Save a plot to a PDF file using `cairo_pdf`
#'
#' @param fn Filename of plot
#' @param p Plot object or function to create plot
#' @param w Width in inches
#' @param h Height in inches
#' @param seed Integer to seed RNG for consistent jittering.
#' @param fun_args List containing arguments to `p` if it's a function and has
#'     required arguments.
#' @param ... Other arguments to `cairo_pdf`
#'
#' @return `NULL`
#' @export
#'
save_plot <- function(fn, p, w, h, seed = NULL, fun_args = NULL, ...) {
    ext <- tail(strsplit(fn, "\\.")[[1]], 1)
    if (ext != "pdf") stop("ERROR: file name extension must be \".pdf\"")
    fn_dir <- dirname(fn)
    if (!dir.exists(fn_dir)) stop("ERROR: `", fn_dir, "` doesn't exist")
    if (!is.null(seed)) set.seed(seed)
    cairo_pdf(filename = fn, width = w, height = h, ...)
    if (is.function(p)) {
        if (is.null(fun_args)) {
            p()
        } else {
            do.call(p, fun_args)
        }
    } else {
        plot(p)
    }
    dev.off()
    invisible(NULL)
}
