#' @title Range of variation
#'
#' @description Selecting the range of differences between the maximum and minimum variation specimen
#'
#' @param procrustes Procrustes data of class \code{"gpagen"}.
#' @param type Which type of coordinates to calculate (see \code{\link{coordinates.differences}} - default is \code{"sperical"}).
#' @param angle Which type of angle to calculate (see \code{\link{coordinates.differences}} - default is \code{"degree"}).
#' @param what Which element from the \code{\link{coordinates.differences}} to use (default is \code{"radius"}).
#' @param ordination Optional, either \code{TRUE} to perform an ordination or directly an ordinated PCA matrix (\code{"prcomp"}) to calculate the range from there.
#' @param axis Optional, if an ordinated matrix is used, which axis (axes) to use. If left empty, all the axes will be used.
#' @param return.ID \code{logical}, whether to return the ID of the max/min specimens or not.
#' @param CI Optional, a value of confidence interval to use (rather than the max/min).
#' 
#' @examples
#' ## Loading the geomorph dataset
#' require(geomorph)
#' data(plethodon)
#' 
#' ## Performing the Procrustes superimposition
#' proc_super <- gpagen(plethodon$land, print.progress = FALSE)
#' 
#' 
#' ## Getting the two most different specimen based on their landmark change radii
#' spec_range <- variation.range(proc_super, return.ID = TRUE)
#' 
#' ## The minimum and maximum specimen
#' spec_range$min.max
#' 
#' ## The range of variation per landmark
#' spec_range$range
#' 
#' 
#' ## Getting the two most different specimen based on the first axis of the ordination
#' variation.range(proc_super, ordination = TRUE, axis = 1, type = "vector", what = "length")
#' 
#' 
#' ## Getting the range variation between specimen using a 95 confidence interval range
#' spec_range95 <- variation.range(proc_super, CI = 0.95, return.ID = TRUE)
#' 
#' ## The absolute maximum and minimum specimens
#' spec_range$min.max
#' 
#' ## The lower and upper 95% range CI specimens
#' spec_range95$min.max
#'
#' @seealso \code{\link{coordinates.difference}}, \code{\link{area.diff}}, \code{\link{rand.test}}
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom stats prcomp

variation.range <- function(procrustes, type = "spherical", angle = "degree", what = "radius", ordination, axis, return.ID = FALSE, CI) {

    ## Sanitizing
    ## procrustes
    check.class(procrustes, "gpagen")

    ## type and angle are dealt with by coordinates.data

    ## what is dealt with by coordinates.area

    ## CI
    if(missing(CI)) {
        do_CI <- FALSE
    } else {
        do_CI <- TRUE
        check.length(CI, 1, msg = " must be on confidence interval in probability.")
        if(CI < 0 || CI > 1) {
            stop("CI must be a percentage between 0 and 1.")
        }
    }

    ## ordination
    if(!missing(ordination)) {

        if(class(ordination) == "logical") {
            if(ordination) {
                do_ordinate <- TRUE

                ## Convert the array in 2D
                array_2d <- two.d.array(procrustes$coords)

                ## measure the tolerance
                tol <- stats::prcomp(array_2d)$sdev^2
                tolerance <- cumsum(tol)/sum(tol)
                tolerance <- length(which(tolerance < 1)) 
                if(length(tolerance) < length(tol)){
                    tolerance <- tolerance + 1
                }
                tolerance <- max(c(tol[tolerance]/tol[1],0.005))

                ## Ordinating the data
                ordination <- stats::prcomp(array_2d, center = TRUE, scale. = FALSE, retx = TRUE, tol = tolerance)

            } else {
                do_ordinate <- FALSE
            }


        } else {
            ## Ordination is a prcomp object
            check.class(ordination, "prcomp")
            do_ordinate <- TRUE
        }

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

        ## Finding the max specimen
        if(do_CI) {
            ## Take the the specimen in the upper CI
            max_specimen <- which(areas == max(areas[which(areas <= quantile(areas, probs = CI))]))  #add abs(areas)?
        } else {
            ## Take the actual max specimen
            max_specimen <- which(areas == max(areas))  #add abs(areas)?
        }

        ## Get the distances from the maximum
        diff_from_max <- coordinates.difference(procrustes$coords[, , -max_specimen], procrustes$coords[, , max_specimen], type = type, angle = angle)

        ## Getting all the areas
        areas_max <- unlist(lapply(diff_from_max, coordinates.area, what = what))

        ## Finding the min specimen
        if(do_CI) {
            ## Take the the specimen in the upper CI
            min_specimen <- which(areas_max == max(areas_max[which(areas_max <= quantile(areas_max, probs = CI))]))  #add abs(areas)?
        } else {
            ## Take the actual max specimen
            min_specimen <- which(areas_max == max(areas_max))  #add abs(areas)?
        }

        ## Get the variation range
        variation_range <- coordinates.difference(procrustes$coords[, , min_specimen], procrustes$coords[, , max_specimen], type = type, angle = angle)[[1]]
    

    } else {

        ## Internal function from geomorph:plotTangentSpace
        get.pc.min.max <- function(axis, what, PCA, GPA, CI) {
            if(length(axis) == 1) {
               output <- arrayspecs(as.vector(t(GPA$consensus)) + c(what(PCA$x[,axis], CI = CI), rep(0, ncol(PCA$x)-length(axis))) %*% t(PCA$rotation), dim(GPA$consensus)[1], dim(GPA$consensus)[2])
            } else {
               output <- arrayspecs(as.vector(t(GPA$consensus)) + c(apply(PCA$x, 2, what, CI = CI), rep(0, ncol(PCA$x)-length(axis))) %*% t(PCA$rotation), dim(GPA$consensus)[1], dim(GPA$consensus)[2])
            }
        }

        ## Get the selector function
        if(do_CI) {
            fun_max <- function(x, CI) return(max(x[which(x <= quantile(x, probs = CI))]))
            fun_min <- function(x, CI) return(min(x[which(x >= quantile(x, probs = 1-CI))]))
        } else {
            CI <- NULL
            fun_max <- function(x, CI) return(max(x))
            fun_min <- function(x, CI) return(min(x))
        }

        ## Applying the method the an ordination
        max_coordinates <- get.pc.min.max(axis = axis, what = fun_max, PCA = ordination, GPA = procrustes, CI = CI)
        dimensions <- dim(max_coordinates)
        max_coordinates <- matrix(max_coordinates, dimensions[1], dimensions[2])
        min_coordinates <- get.pc.min.max(axis = axis, what = fun_min, PCA = ordination, GPA = procrustes, CI = CI)
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