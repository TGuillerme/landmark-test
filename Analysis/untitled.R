#' @title Area difference
#'
#' @description Measure the absolute area difference between two distribution
#'
#' @param 
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export


area.diff <- function(x, y) {
    ## Sort both distributions (y axis)
    y_x <- sort(x)
    y_y <- sort(y)

    ## Get the x coordinates
    x <- 1:max(length(x), length(y))
    
    ## Calculate the x~y area
    x_area <- sum(diff(x_x) * zoo::rollmean(y_x,2))
    y_area <- sum(diff(x_y) * zoo::rollmean(y_y,2))
    
    return(x_area - y_area)
}
