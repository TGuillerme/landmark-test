#' @title Area difference
#'
#' @description Measure the area difference between two ranked distribution
#'
#' @param x,y the two distributions to compare.
#' @param rarefy Optional, if \code{x} is \neq \code{y}, how many rarefaction samples to use (\code{defaut = 500})
#' @param cent.tend Optional, if \code{x} is \neq \code{y}, which central tendency to use (\code{defaut = mean})
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

area.diff <- function(x, y, rarefy = 500, cent.tend = mean) {
    ## Sort both distributions (y axis)
    y_x <- sort(x)
    y_y <- sort(y)

    ## Getting the vectors lengths
    length_x <- length(x)
    length_y <- length(y)

    if(length_x == length_y) {
        do_rarefy <- FALSE
        ## Samples have the same size
        x_axis <- 1:length_x
    } else {
        do_rarefy <- TRUE
        ## Samples have different sizes
        if(length_x > length_y) {
            x_rare <- TRUE
            ## x > y
            x_rarefied <- lapply(replicate(rarefy, sample(x, length_y), simplify = FALSE), sort)
            y_rarefied <- y_y
            x_axis <- 1:length_y
        } else {
            x_rare <- FALSE
            ## x < y
            y_rarefied <- lapply(replicate(rarefy, sample(y, length_x), simplify = FALSE), sort)
            x_rarefied <- y_x
            x_axis <- 1:length_x
        }
    }

    if(!do_rarefy) {
        ## Calculate the x~y area
        x_area <- sum(diff(x_axis) * zoo::rollmean(y_x,2))
        y_area <- sum(diff(x_axis) * zoo::rollmean(y_y,2))
        
    } else {

        if(x_rare) {
            x_area <- cent.tend(unlist(lapply(x_rarefied, function(X, x_axis) return(sum(diff(x_axis) * zoo::rollmean(X,2))), x_axis)))
            y_area <- sum(diff(x_axis) * zoo::rollmean(y_y,2))
        } else {
            x_area <- sum(diff(x_axis) * zoo::rollmean(y_x,2))
            y_area <- cent.tend(unlist(lapply(y_rarefied, function(X, x_axis) return(sum(diff(x_axis) * zoo::rollmean(X,2))), x_axis)))
        }

    }

    return(x_area - y_area)
}
