#' @title Coordinates area
#'
#' @description Calculates the area of the coordinates differences
#'
#' @param coordinates.diff The \code{matrix} of coordinate differences.
#' @param what Which element of the coordinate differences to use (can be the \code{numeric} value of the column or the name as \code{character}).
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom zoo rollmean


coordinates.area <- function(coordinate_diff, what = "radius") {
    y <- sort(coordinate_diff[, what])
    x <- 1:length(y)
    return(sum(diff(x) * zoo::rollmean(y,2)))
}

##### TESTS
test <- FALSE
if(test){

    context("coordinates.area")

}