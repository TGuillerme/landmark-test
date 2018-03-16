#' @title Add the rarefaction plot
#'
#' @description Add the rarefaction plot to a called \code{\link[ade4]{plot.randtest}} plot
#'
#' @param x The \code{randtest} object
#' @param ... Any optional parameters to be passed to \code{\link[ade4]{plot.randtest}} or \code{\link[graphics]{plot}}.
#' 
#' @examples
#' ## Loading the geomorph dataset
#' require(geomorph)
#' data(plethodon)
#' 
#' ## Performing the Procrustes superimposition
#' proc_super <- gpagen(plethodon$land, print.progress = FALSE)
#' 
#' ## Getting the two most different specimen based on their landmark change radii
#' var_range <- variation.range(proc_super)
#' 
#' \dontrun{
#' ## Rarefying the area difference to 4 elements without testing the parameter
#' rarefy_test <- rand.test(var_range[, "radius"], random_part, rarefaction = 5,
#'                          test = area.diff)
#' 
#' ## Plotting the results (central tendency)
#' plot(rarefy_test)
#' 
#' ## Add the rarefied observations to the plot
#' add.rare.plot(rarefy_test)
#' }
#' 
#' @seealso \code{\link{rand.test}}
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom graphics hist lines points

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
    graphics::hist(rarefied_distribution, col = dots$col, add = TRUE)

    ## Add the central tendency as plot.randtest (code directly from: https://github.com/cran/ade4/blob/master/R/randtest.R)
    h0 <- x$plot$hist
    y0 <- max(h0$counts)
    graphics::lines(c(x$obs, x$obs), c(y0/2, 0))
    graphics::points(x$obs, y0/2, pch = 18, cex = 2)

    return(invisible())
}