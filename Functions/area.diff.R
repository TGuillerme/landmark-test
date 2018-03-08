#' @title Area difference
#'
#' @description Measure the area difference between two ranked distribution
#'
#' @param x,y the two distributions to compare.
#' @param rarefy Optional, if \code{x} is \neq \code{y}, how many rarefaction samples to use. If left empty a default number of replicates between 100 and 1000 is used (see details).
#' @param cent.tend Optional, if \code{x} is \neq \code{y}, which central tendency to use (\code{defaut = mean})
#' 
#' @details
#' The number of replicates is chosen based on the variance of the distribution to rarefy using the Silverman's rule of thumb (Silverman 1986, pp.48, eqn 3.31) for choosing the bandwidth of a Gaussian kernel density estimator multiplied by 1000 with a result comprised between 100 and 1000.
#' For example, for rarefying \code{x}, \code{rarefy = round(stats::bw.nrd0(x) * 1000)}.
#' With \code{100 \leq rarefy \leq 1000}.
#' 
#' @references
#' Silverman, B. W. (1986) Density Estimation. London: Chapman and Hall.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom zoo rollmean
#' @importFrom stats bw.nrd0

area.diff <- function(x, y, rarefy, cent.tend = mean) {
    ## Sort both distributions (y axis)
    y_x <- sort(x, decreasing = TRUE)
    y_y <- sort(y, decreasing = TRUE)

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

            ## If rarefying is missing sample 1000 times the sd with a max of 500 and a min of 100
            if(missing(rarefy)) {
                rarefy <- round(stats::bw.nrd0(x) * 1000)
                rarefy <- ifelse(rarefy > 1000, 1000, rarefy)
                rarefy <- ifelse(rarefy < 100, 100, rarefy)
            }

            ## x > y
            x_rarefied <- lapply(replicate(rarefy, sample(x, length_y), simplify = FALSE), sort)
            y_rarefied <- y_y
            x_axis <- 1:length_y
        } else {
            x_rare <- FALSE

            ## If rarefying is missing sample 1000 times the sd with a max of 500 and a min of 100
            if(missing(rarefy)) {
                rarefy <- round(stats::bw.nrd0(y) * 1000)
                rarefy <- ifelse(rarefy > 1000, 1000, rarefy)
                rarefy <- ifelse(rarefy < 100, 100, rarefy)
            }

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
