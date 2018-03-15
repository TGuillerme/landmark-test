#' @title Coordinates differences
#'
#' @description Calculates coordinate differences (e.g. from procrustes superimpositions)
#'
#' @param coordinates An \code{array}, \code{list} or \code{matrix} of coordinates to be compared to the reference.
#' @param reference A \code{matrix} of reference coordinates. If missing and \code{coordinates} is an \code{array} or \code{list}, the first element of \code{coordinates} is used.
#' @param type the type of coordinates to output: can be \code{"cartesian"} (x0,y0,x1,y1,... format), or  \code{"spherical"} (radius, polar, azimuth).
# \code{"vector"} (length, angle(s))
#' @param angle optional, whether display angles in radian (\code{"radian"} - default) or in degrees (\code{"degree"}).
#' 
#' @examples
#' ## Loading the geomorph dataset
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
#' @seealso \code{\link{variation.range}}
#' 
#' @author Thomas Guillerme
#' @export

coordinates.difference <- function(coordinates, reference, type = "cartesian", angle = "radian") {

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
            coordinates <- sapply(1:dim(coordinates)[3], function(x, coordinates) return(list(coordinates[,,x])), coordinates)
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

    ## Method
    type <- tolower(type)
    check.method(type, c("cartesian", "spherical"), msg = "type") #"vector"

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


        # fun.angle <- function(one_row, axis, dimension) {
        #     return(
        #         acos( ( sqrt((one_row[-c(1:dimension)][axis] - one_row[1:dimension][axis])^2) ) / sqrt(sum((one_row[-c(1:dimension)]-one_row[1:dimension])^2)) )
        #         )
        # }

        # library(testthat)


        # expect_values <- c()


        # expect_equal(
        # fun.angle(c(1,1,2,2), 1, 2)*180/pi
        # , 45)
        # expect_equal(
        # fun.angle(c(1,1,-2,-2), 1, 2)*180/pi
        # , 45)
        # expect_equal(
        # fun.angle(c(1,1,2,1), 1, 2)*180/pi
        # , 0)
        # expect_equal(
        # fun.angle(c(1,1,2,1), 2, 2)*180/pi
        # , 90)    
        # expect_equal(
        # fun.angle(c(1,1,1,2), 1, 2)*180/pi
        # , 90)
        # expect_equal(
        # fun.angle(c(1,1,1,2), 2, 2)*180/pi
        # , 0)

        ## Calculate the angles
        output <- apply(one_coordinate, 1, fun.angle, axis = axis, dimension = dimension)
        
        ## Convert degrees
        if(degree) { output <- output*180/pi}




        return(output)
    }


    ## Vector coordinates
    # if(type == "vector") {
    #     ## length
    #     length <- lapply(coordinates, euclidean.distance, dimensions[2])

    #     ## angle
    #     if(dimensions[2] == 2) {
    #         ## Get the x angle
    #         angle <- lapply(coordinates, get.angle, axis = 1, dimensions[2], degree = degree)

    #         ## Direction
    #         direction <- lapply(angle, function(x) return(cos(x)/sin(x)))

    #     } else {
    #         ## Get the x and y angles
    #         angle_x <- lapply(coordinates, get.angle, axis = 1, dimensions[2], degree = degree)
    #         angle_y <- lapply(coordinates, get.angle, axis = 2, dimensions[2], degree = degree)
    #         ## Combine both
    #         angle <- mapply(cbind, angle_x, angle_y, SIMPLIFY = FALSE)

    #         ## Direction
    #         direction <- lapply(angle_x, function(x) return(cos(x)/sin(x)))

    #     }

    #     ## Combine the coordinates
    #     coordinates <- mapply(cbind, length, angle, SIMPLIFY = FALSE)
    #     coordinates <- mapply(cbind, coordinates, direction, SIMPLIFY = FALSE)

    #     if(dimensions[2] == 2) {
    #         coordinates <- lapply(coordinates, function(x) {colnames(x) <- c("length", "angle", "direction") ; return(x)})
    #     } else {
    #         coordinates <- lapply(coordinates, function(x) {colnames(x) <- c("length", "angle_x", "angle_y", "direction") ; return(x)})
    #     }

    # }

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