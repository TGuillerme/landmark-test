#' @title Random test
#'
#' @description Performs a random test
#'
#' @param distribution A \code{numeric} distribution to draw from.
#' @param subset A \code{numeric} vector of the elements of the distribution to test.
#' @param test The test to apply (a \code{function}).
#' @param replicates A \code{numeric} value for the number of replicates.
#' @param rarefaction A \code{numeric} value for rarefying the subset.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

rand.test <- function(distribution, subset, test, replicates, rarefaction) {

    return()
}


# ## Randomly sample
# rand.sample <- function(max_min, selected_lan) {
#     ## Randomly sampling
#     random_select <- sample(1:length(max_min), length(selected_lan))

#     ## Wilcox test
#     wilcox <- wilcox.test(max_min[random_select], max_min[-random_select])

#     ## Bhattaacharrya
#     bhatt <- bhatt.coeff(max_min[random_select], max_min[-random_select])

#     ## Results
#     return(c("BC" = bhatt, wilcox$statistic, wilcox$p.value))
# }

# ## Doing it on 100 bootstraps
# test <- replicate(100, rand.sample(first_axis_max_min[[1]][,1], selected_land))


##### TESTS
test <- FALSE
if(test){
    context("rand.test")

    #Test
    test_that("rand.test works", {

        ## Output

    })
}