#' @title Range of variation
#'
#' @description Selecting the range of differences between the maximum and minimum variation specimen
#'
#' @param procrustes Procrustes data of class \code{"gpagen"}.
#' @param type Which type of coordinates to calculate (see \code{\link{coordinates.differences}} - default is \code{"spherical"}).
#' @param angle Which type of angle to calculate (see \code{\link{coordinates.differences}} - default is \code{"degree"}).
#' @param what Which element from the \code{\link{coordinates.differences}} to use (default is \code{"radius"}).
#' @param ordination Optional, an ordinated PCA matrix (\code{"prcomp"}) to calculate the range from there.
#' @param axis Optional, if an ordinated matrix is used, which axis (axes) to use. If left empty, all the axes will be used.
#' @param return.ID \code{logical}, whether to return the ID of the max/min specimens or not.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

variation.range <- function(procrustes, type = "spherical", angle = "degree", what = "radius", ordination, axis, return.ID = FALSE) {

    ## Sanitizing
    ## procrustes
    check.class(procrustes, "gpagen")

    ## type and angle are dealt with by coordinates.data

    ## what is dealt with by coordinates.area

    ## ordination
    if(!missing(ordination)) {
        check.class(ordination, "prcomp")
        do_ordinate <- TRUE
    } else {
        do_ordinate <- FALSE
    }

    ## axis
    if(do_ordinate) {
        if(missing(axis)) {
            axis <- 1:ncol(ordination$x)
        } else {
            check.class(axis, "numeric")
            if(min(axis) < 1 || max(axis) > ncol(ordination$x)) {
                stop("Wrong axis number.")
            }
        }
    }
    
    ## return.ID
    check.class(return.ID, "logical")

    ## Applying the method to the Procrustes
    if(!do_ordinate) {

        ## Get the distances from the consensus
        diff_consensus <- coordinates.difference(procrustes$coords, procrustes$consensus, type = type, angle = angle)

        ## Get the volume of change for each element (area under the curve)
        areas <- unlist(lapply(diff_consensus, coordinates.area, what = what))

        ## Finding the specimen the most different from the consensus
        max_specimen <- which(areas == max(areas))

        ## Get the distances from the maximum
        diff_from_max <- coordinates.difference(procrustes$coords, procrustes$coords[,,max_specimen], type = type, angle = angle)

        ## Getting all the areas
        areas_max <- unlist(lapply(diff_from_max, coordinates.area, what = what))

        ## Finding the specimen the most different from the consensus
        min_specimen <- which(areas_max == max(areas_max))

        ## Get the variation range
        variation_range <- coordinates.difference(procrustes$coords[,,min_specimen], procrustes$coords[,,max_specimen], type = type, angle = angle)[[1]]
    

    } else {

        ## Internal function from geomorph:plotTangentSpace
        get.pc.min.max <- function(axis, what, PCA, GPA) {
            if(length(axis) == 1) {
               output <- geomorph::arrayspecs(as.vector(t(GPA$consensus)) + c(what(PCA$x[,axis]), rep(0, ncol(PCA$x)-length(axis))) %*% t(PCA$rotation), dim(GPA$consensus)[1], dim(GPA$consensus)[2])
            } else {
               output <- geomorph::arrayspecs(as.vector(t(GPA$consensus)) + c(apply(PCA$x, 2, what), rep(0, ncol(PCA$x)-length(axis))) %*% t(PCA$rotation), dim(GPA$consensus)[1], dim(GPA$consensus)[2])
            }
        }

        ## Applying the method the an ordination
        max_coordinates <- get.pc.min.max(axis = axis, what = max, PCA = ordination, GPA = procrustes)
        dimensions <- dim(max_coordinates)
        max_coordinates <- matrix(max_coordinates, dimensions[1], dimensions[2])
        min_coordinates <- get.pc.min.max(axis = axis, what = min, PCA = ordination, GPA = procrustes)
        min_coordinates <- matrix(min_coordinates, dimensions[1], dimensions[2])

        ## Get the variation range
        variation_range <- coordinates.difference(min_coordinates, max_coordinates, type = type, angle = angle)[[1]]
    }
    
    if(!do_ordinate && return.ID){   
        return(list("range" = variation_range, "min.max" = c(min_specimen, max_specimen)))
    } else {
        return(variation_range)
    }
}