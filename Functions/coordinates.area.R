#' @title Coordinates area
#'
#' @description Calculates the area of the coordinates differences
#'
#' @param data The \code{matrix} of coordinate differences.
#' @param what Which element of the coordinate differences to use (can be the \code{numeric} value of the column or the name as \code{character}).
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom zoo rollmean


coordinates.area <- function(data, what = "radius") {

    ## Sanitising
    class_what <- class(what)
    if(class_what %in% c("numeric", "integer")) {
        if(what > ncol(data)) {
            stop(paste0("what argument cannot be greater than the number of columns available in data (", ncol(data) ,")."))
        }
    } else {
        if(class_what == "character") {
            if(!(what %in% colnames(data))) {
                stop(paste0(what, " not found in data column names."))
            }
        } else {
            stop("what argument must be either of class 'character' or 'numeric'.")
        }
    }

    ## Get the y coordinates
    y <- sort(data[, what])

    ## Get the x coordinates
    x <- 1:length(y)
    
    ## Calculate the x~y area
    return(sum(diff(x) * zoo::rollmean(y,2)))
}

# ##### TESTS
# test <- FALSE
# if(test){

#     context("coordinates.area")

# }