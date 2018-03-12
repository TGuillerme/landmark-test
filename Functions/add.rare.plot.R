#' @title Add the rarefaction plot
#'
#' @description Add the rarefaction plot to a called randtest plot
#'
#' @param x The \code{randtest} object
#' @param ... Any optional parameters to be passed to \code{plot.randtest} or \code{plot}
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

add.rare.plot <- function(x, ...) {

    ## Sanitizing
    check.class(x, "randtest")
    if(is.null(x$observed)) {
        stop("x must have an $observed element returned from rand.test().")
    } else {
        if(length(x$observed) == 1) {
            stop("x must have an $observed element returned from rand.test().")
        } else {
            rarefied_distribution <- unlist(x$observed)
        }
    }

    ## optionals
    dots <- list(...)
    if(is.null(dots$col)) {
        dots$col <- "white"
    }

    ## Add the histogram
    hist(rarefied_distribution, col = dots$col, add = TRUE)

    ## Add the central tendency as plot.randtest (code directly from: https://github.com/cran/ade4/blob/master/R/randtest.R)
    h0 <- x$plot$hist
    y0 <- max(h0$counts)
    lines(c(x$obs, x$obs), c(y0/2, 0))
    points(x$obs, y0/2, pch = 18, cex = 2)

    return(invisible())
}