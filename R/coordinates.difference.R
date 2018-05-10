#' @title Coordinates differences
#'
#' @description Calculates coordinate differences (e.g. from procrustes superimpositions)
#'
#' @param coordinates An \code{array}, \code{list} or \code{matrix} of coordinates to be compared to the reference.
#' @param reference A \code{matrix} of reference coordinates. If missing and \code{coordinates} is an \code{array} or \code{list}, the first element of \code{coordinates} is used.
#' @param type the type of coordinates to output: can be \code{"cartesian"} (x0,y0,x1,y1,... format), \code{"spherical"} (radius, polar, azimuth) or \code{"vector"} (length, angle(s)).
#' @param angle optional, whether display angles in radian (\code{"radian"} - default) or in degrees (\code{"degree"}).
#' @param absolute.distance \code{logical}, when using \code{"vector"}, whether to use the absolute distance (from the centroid) or the relative one (from the reference landmark). See details.
#' 
#' @details
#' When using \code{type = "vector"} with \code{absolute.distance = TRUE}, the distance between two landmarks A and A' is calculated as d(0,A') - d(0,A) where 0 is the centroid of the shape to analysis.
#' A positive absolute distance means that A' is further away from the centroid than A.
#' A negative absolute distance means that A' closer to the centroid than A.
#' When using \code{absolute.distance = FALSE}, the distance is calculated as d(A,A') and is always positive.
#' 
#' @examples
#' ## Loading the geomorph dataset
#' require(geomorph)
#' data(plethodon)
#' 
#' ## Performing the Procrustes superimposition
#' proc_super <- geomorph::gpagen(plethodon$land, print.progress = FALSE)
#' 
#' ## Getting the coordinates differences from the consensus
#' cartesian_diff <- coordinates.difference(proc_super$coords, proc_super$consensus)
#' 
#' ## The coordinates of the differences between the first specimen and the consensus
#' head(cartesian_diff[[1]])
#' 
#' ## Getting the spherical coordinates difference between the two first specimen
#' coordinates.difference(proc_super$coords[, , 1], proc_super$coords[, , 2],
#'                        type = "spherical", angle = "degree")
#' 
#' ## Getting the vector coordinates for the same specimen in relative distance
#' coordinates.difference(proc_super$coords[, , 1], proc_super$coords[, , 2],
#'                        type = "vector", angle = "degree", absolute.distance = FALSE)
#' 
#' ## Getting the vector same coordinates in absolute distances
#' coordinates.difference(proc_super$coords[, , 1], proc_super$coords[, , 2],
#'                        type = "vector", angle = "degree")
#' @seealso \code{\link{variation.range}}
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom geometry dot

coordinates.difference <- function(coordinates, reference, type = "cartesian", angle = "radian", absolute.distance = TRUE) {

#coordinates.difference(procrustes$coords, procrustes$coords[,,max_specimen], type = type, angle = angle)

    ## Sanitizing

    ## Coordinates
    coordinates_class <- check.class(coordinates, c("array", "list", "matrix"))
    if(coordinates_class == "matrix") {
        ## Get dimensions
        dimensions <- dim(coordinates)
        coordinates <- list(coordinates)
    } else {
        if(coordinates_class == "array") {
            ## Convert array into list
            if(!is.null(attributes(coordinates)$dimnames[[3]])) {
                coordi_names <- attributes(coordinates)$dimnames[[3]]
            } else {
                coordi_names <- 1:dim(coordinates)[3]
            }
            coordinates <- sapply(1:dim(coordinates)[3], function(x, coordinates) return(list(coordinates[,,x])), coordinates)
            names(coordinates) <- coordi_names
        }
        ## Get dimensions
        dimensions <- unique(lapply(coordinates, dim))
        if(length(dimensions) != 1) {
            stop("Elements in coordinates don't have the same dimensions!")
        } else {
            dimensions <- unlist(dimensions)
        }
    }

    ## Reference
    if(missing(reference) && coordinates_class != "matrix") {
        ## Default reference
        reference <- coordinates[[1]]
    } else {
        check.class(reference, "matrix")
        ## Check the dimensions
        ref_dimensions <- dim(reference)
        if(!all(ifelse(ref_dimensions == dimensions, TRUE, FALSE))) stop("reference has not the same dimensions as coordinates!")
    }

    ## Check if the coordinates have names
    add_names <- ifelse(is.null(names(coordinates)), FALSE, TRUE)
    if(add_names) {
        coordi_names <- names(coordinates)
    } else {
        ## Check if the names are not in an array dimnames attributes
        if(!is.null(attributes(coordinates)$dimnames[[3]])) {
            add_names <- TRUE
            coordi_names <- attributes(coordinates)$dimnames[[3]]
        }
    }

    ## Method
    type <- tolower(type)
    check.method(type, c("cartesian", "spherical", "vector"), msg = "type")

    ## Angle
    angle <- tolower(angle)
    check.method(angle, c("radian", "degree"), msg = "angle")
    degree <- ifelse(angle == "degree", TRUE, FALSE)

    ## Getting the vector coordinates
    get.coord <- function(one_coordinate, reference, dimension) {
        ## Get the coordinates
        coordinate_matrix <- cbind(reference, one_coordinate)
        ## Name them
        if(dimensions[2] < 27) {
            name_avail <- c(letters[24:26],letters[23:1])
        } else {
            name_avail <- paste0("d",seq(1:dimensions[2]),"_")
        }
        colnames(coordinate_matrix) <- c(paste0(name_avail[1:dimension], 0), paste0(name_avail[1:dimension], 1))
        return(coordinate_matrix)
    }

    ## Calculate the coordinates
    coordinates <- lapply(coordinates, get.coord, reference, dimension = dimensions[2])

    ## Function for getting different coordinates
    euclidean.distance <- function(one_coordinate, dimension) {
        fun.dist <- function(one_row, dimension) {
            return(sqrt(sum((one_row[-c(1:dimension)]-one_row[1:dimension])^2)))
        }
        return(apply(one_coordinate, 1, fun.dist, dimension))
    }

    get.angle <- function(one_coordinate, axis, dimension, degree) {

        ## Absolute version
        fun.angle <- function(one_row, axis, dimension) {
            return(
                acos( ( sqrt((one_row[-c(1:dimension)][axis] - one_row[1:dimension][axis])^2) ) / sqrt(sum((one_row[-c(1:dimension)]-one_row[1:dimension])^2)) )
                )
        }

        ## Calculate the angles
        output <- apply(one_coordinate, 1, fun.angle, axis = axis, dimension = dimension)
        
        ## Convert degrees
        if(degree) { output <- output*180/pi}

        return(output)
    }

    get.vector.diffs <- function(one_coordinate, dimension, angle, absolute.distance = TRUE) {

        ## Transform coordinates into a vector
        coord.to.vector <- function(coords, dimension) {
            mapply(function(a, b) return(a - b), coords[1:dimension], coords[-c(1:dimension)])
        }

        ## Get centroid
        centroid <- apply(one_coordinate[,1:dimension], 2, mean)

        ## Get the vectors for the reference
        reference <- cbind(matrix(rep(centroid, nrow(one_coordinate)), ncol = dimension, byrow = TRUE), one_coordinate[,1:dimension])
        reference <- apply(reference, 1, coord.to.vector, dimension)
        #rownames(reference) <- paste0("u", 1:dimension)

        ## Get the vector for the observed
        observed <- cbind(matrix(rep(centroid, nrow(one_coordinate)), ncol = dimension, byrow = TRUE), one_coordinate[,-c(1:dimension)])
        observed <- apply(observed, 1, coord.to.vector, dimension)
        #rownames(observed) <- paste0("v", 1:dimension)

        ## Get the vector lengths
        ref_length <- apply(reference, 2, function(X) return(sqrt(sum(X^2))))
        obs_length <- apply(observed, 2, function(X) return(sqrt(sum(X^2))))

        ## Dot product of the vectors
        dot_prod <- geometry::dot(reference, observed)

        ## Get the angle
        angles <- acos(dot_prod / (ref_length * obs_length))

        if(angle == "degree") {
            angles <- angles * 180 / pi
        }

        ## Get the length difference
        if(absolute.distance == TRUE) {
            length <- obs_length - ref_length
        } else {
            length <- euclidean.distance(one_coordinate, dimension)
        }

        return(matrix(c(length, angles), ncol = 2, byrow = FALSE, dimnames = list(c(),c("length", "angle"))))
    }

    ## Spherical coordinates
    if(type == "spherical") {

        ## Get the length (radial)
        radius <- lapply(coordinates, euclidean.distance, dimensions[2])

        ## Get the azimuth
        azimuth <- lapply(coordinates, get.angle, axis = 1, dimensions[2], degree = degree)

        if(dimensions[2] != 2) {
            ## Get the polar
            polar <- lapply(coordinates, get.angle, axis = 3, dimensions[2], degree = degree)
        }

        ## Combine the coordinates
        coordinates <- mapply(cbind, radius, azimuth, SIMPLIFY = FALSE)

        if(dimensions[2] == 2) {
            coordinates <- lapply(coordinates, function(x) {colnames(x) <- c("radius", "azimuth") ; return(x)})
        } else {
            coordinates <- mapply(cbind, coordinates, polar, SIMPLIFY = FALSE)
            coordinates <- lapply(coordinates, function(x) {colnames(x) <- c("radius", "azimuth", "polar") ; return(x)})
        }
    }

    ## Vector coordinates
    if(type == "vector") {

        ## Get the vector coordinates for each coordinates
        coordinates <- lapply(coordinates, get.vector.diffs, dimension = dimensions[2], angle = angle, absolute.distance = absolute.distance)

        # test <- list()
        # for(coord in 1:length(coordinates)) {
        #     cat(paste0(coord, "\n"))
        #     test[[coord]] <- get.vector.diffs(coordinates[[coord]], dimension = dimensions[2], angle = angle, absolute.distance = absolute.distance)
        # }
    }

    if(add_names) {
        names(coordinates) <- coordi_names
    }

    return(coordinates)
}


# ##### TESTS
# test <- FALSE
# if(test){
#     #TESTING bhatt.coeff

#     context("coordinates.difference")

#     #Test
#     test_that("correct output", {

#         ## Make data set
#         data(plethodon)
#         test_output <- geomorph::gpagen(plethodon$land, PrinAxes = FALSE, print.progress = FALSE)
#         test_out <- coordinates.difference()

#     })
# }


# stop("DEBUG in coordinates difference")

# ## Get the ranges
# range_1 <- apply(proc_super$coords[,, 33], 2, range)
# range_2 <- apply(proc_super$coords[,, 14], 2, range)
# xlim <- ylim <- range(as.vector(range_1), as.vector(range_2))

# par(bty = "n")
# plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "")

# ## Adding the consensus shape
# points(proc_super$consensus, pch = 19, col = "grey")
# centroid <- apply(proc_super$consensus, 2, mean)
# points(x = centroid[1], y = centroid[2], pch = 13, cex = 2)
# for(land in 1:nrow(proc_super$consensus)) {
#     lines(x = c(centroid[1], proc_super$consensus[land,1]), y = c(centroid[2], proc_super$consensus[land,2]), col = "grey", lty = 3)
# }

# ## Adding the min var
# point_val <- 33
# col_val <- "blue"
# points(proc_super$coords[,, point_val], pch = 19, col = col_val)
# for(land in 1:nrow(proc_super$consensus)) {
#     ## Adding the distance from centroid
#     lines(x = c(centroid[1], proc_super$coords[land,1, point_val]), y = c(centroid[2], proc_super$coords[land,2, point_val]), col = col_val, lty = 3)
#     ## Adding the magnitude
#     lines(x = c(proc_super$consensus[land,1], proc_super$coords[land,1, point_val]), y = c(proc_super$consensus[land,2], proc_super$coords[land,2, point_val]), col = col_val, lty = 1)
# }

# ## Adding the max var
# point_val <- 14
# col_val <- "red"
# points(proc_super$coords[,, point_val], pch = 19, col = col_val)
# for(land in 1:nrow(proc_super$consensus)) {
#     ## Adding the distance from centroid
#     lines(x = c(centroid[1], proc_super$coords[land,1, point_val]), y = c(centroid[2], proc_super$coords[land,2, point_val]), col = col_val, lty = 3)
#     ## Adding the magnitude
#     lines(x = c(proc_super$consensus[land,1], proc_super$coords[land,1, point_val]), y = c(proc_super$consensus[land,2], proc_super$coords[land,2, point_val]), col = col_val, lty = 1)
# }

# ## Add the difference between max and min
# point_val <- c(33, 14)
# col_val <- "black"
# for(land in 1:nrow(proc_super$consensus)) {
#     ## Adding the distance from centroid
#     lines(x = c(proc_super$coords[land,1, point_val[1]], proc_super$coords[land,1, point_val[2]]), y = c(proc_super$coords[land,2, point_val[1]], proc_super$coords[land,2, point_val[2]]), col = col_val, lty = 1)

# }


