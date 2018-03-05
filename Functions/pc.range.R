#' @title PC range
#'
#' @description Gets the coordinates of the PC max and min values
#'
#' @param axis One or more axes to get the coordinates from (\code{numeric}).
#' @param what Which element to get from the PC axis (\code{max}, \code{min}, \code{mean}, etc.)
#' @param PCA An ordinated matrix (must be ordinated using Principal Coordinates Analysis).
#' @param GPA The coordinates of the Procrustes superimposition.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom geomorph arrayspecs

## This is what ploTangentSpace does individually
pc.range <- function(axis, what, PCA, GPA) {
    if(length(axis) == 1) {
       output <- geomorph::arrayspecs(as.vector(t(GPA$consensus)) + c(what(PCA$pc.scores[,axis]),
                                      rep(0, ncol(PCA$pc.scores)-length(axis))) %*% t(PCA$rotation),
                                      dim(GPA$consensus)[1], dim(GPA$consensus)[2])
    } else {

        ## To check with Emma for this one

       output <- geomorph::arrayspecs(as.vector(t(GPA$consensus)) + c(apply(PCA$pc.scores, 2, what),
                                      rep(0, ncol(PCA$pc.scores)-length(axis))) %*% t(PCA$rotation),
                                      dim(GPA$consensus)[1], dim(GPA$consensus)[2])
    }
}

##### TESTS
test <- FALSE
if(test){
    context("pc.range")
}