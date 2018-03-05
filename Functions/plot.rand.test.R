#' @title Plotting the random test results
#'
#' @description Plotting the random test results
#'
#' @param x The results of \code{rand.test}.
#' @param ... any optional arguments to be passed to \code{plot}.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

plot.rand.test <- function(test, obs, which, xlab = "") {
    densi <- hist(test[which,], plot = FALSE)
    plot(NULL, main = "", xlab = xlab, xlim = range(c(densi$breaks, obs)), ylim = range(c(0, densi$counts)), ylab = "Density")
    hist(test[which,], add = TRUE)
    abline(v = obs, lwd = 2, lty = 2)
}

##### TESTS
test <- FALSE
if(test){

    context("plot.rand.test")

}