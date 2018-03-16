#' @title Random test
#'
#' @description Performs a random test
#'
#' @param distribution A \code{numeric} distribution to draw from.
#' @param subset A \code{numeric} vector of the elements of the distribution to test.
#' @param test The test to apply (a \code{function}).
#' @param replicates A \code{numeric} value for the number of replicates (\code{default = 100}).
#' @param resample \code{logical} wether to resample the full distribution (\code{TRUE}) or the distribution without the subset (\code{FALSE}).
#' @param rarefaction Optional, a \code{numeric} value for rarefying the subset.
#' @param test.parameter Optional, whether to test the calculated parameter (\code{TRUE}) or not (\code{FALSE} - default).
#' @param parameter Optional,  if the parameter is tested, which parameter to select (can be left empty if \code{test} outputs a single value or a named object containing a \code{statistic} element).
#' @param alternative Optional, if the parameter is tested, what is the alternative hypothesis. Can be \code{"two-sided"} (default), \code{"greater"} or \code{"lesser"}.
#' @param ... Any optional arguments to be passed to \code{test}.
#' 
#' @details First, the subset of the distribution is compared to the whole distribution (observed difference).
#' Second, random equally sized subsets are compared to the whole distribution (random differences).
#' If the observed difference falls out of the random differences distribution, the differences are significant.
#' 
#' If \code{test.parameter} is set to \code{TRUE}, the code tests whether the resulting test parameters between the observed subset and the random ones are significantly different (base on the same procedure as in \code{link[ade4]{rantest}}).
#' 
#' @return
#' This function returns a \code{"randtest"} object that can be passed to the generic S3 functions \code{\link[randtest]{print.randtest}} or \code{\link[randtest]{plot.randtest}}.
#' The output also contains to extra elements \code{output$observed} and \code{output$random} containing the raw results of respectively the observed and random tests.
#' 
#' @examples
#' ## Loading the geomorph dataset
#' require(geomorph)
#' data(plethodon)
#' 
#' ## Performing the Procrustes superimposition
#' proc_super <- gpagen(plethodon$land, print.progress = FALSE)
#' 
#' ## Getting the two most different specimen based on their landmark change radii
#' var_range <- variation.range(proc_super)
#' 
#' set.seed(1)
#' 
#' ## Selecting 6 random landmarks
#' random_part <- sample(1:nrow(var_range), 6)
#' 
#' ## Testing the difference between the two sets of landmarks
#' stats::t.test(var_range[random_part, "radius"], var_range[-random_part, "radius"])
#' 
#' ## Testing whether this difference is expected by chance
#' random_test <- rand.test(var_range[, "radius"], random_part, test = stats::t.test,
#'                         test.parameter = TRUE)
#'
#' ## Summarising the results
#' random_test
#' 
#' ## Plotting the results
#' plot(random_test)
#' 
#' ## Rarefying the area difference to 4 elements without testing the parameter
#' rarefy_test <- rand.test(var_range[, "radius"], random_part, rarefaction = 5,
#'                          test = t.test)
#' plot(rarefy_test)
#'
#' @seealso \code{link[ade4]{randtest}}, \code{\link{bootstrap.test}}
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom stats sd var
#' @importFrom graphics hist

rand.test <- function(distribution, subset, test, replicates = 100, resample = TRUE, rarefaction, test.parameter = FALSE, parameter, alternative = "two-sided", ...) {

    match_call <- match.call()

    ## Sanitizing
    ## Distribution and subset
    check.class(distribution, "numeric")
    check.class(subset, c("numeric", "integer"))
    if(length(subset) > length(distribution)) {
        stop("The subset is bigger than the whole distribution.")
    }

    ## Test
    check.class(test, "function")

    ## Replicates
    check.class(replicates, "numeric")
    check.length(replicates, 1, msg = " must be a single numeric value.")
    if(replicates < 1) {
        stop("At least one replicate must be run.")
    }

    ## Resample
    check.class(resample, "logical")

    ## Rarefaction
    if(missing(rarefaction)) {
        is_rarefied <- FALSE
    } else {
        is_rarefied <- TRUE
        check.class(rarefaction, c("numeric", "integer"))
        check.length(rarefaction, 1, msg = " must be a single numeric value.")
        if(length(rarefaction) > length(subset)) {
            stop("The rarefaction values is bigger than the subset size.")
        }
    }

    ## test.parameter
    check.class(test.parameter, "logical")

    ## Testing the parameter (if any)
    test_parameter <- FALSE
    ## Testing the test output
    testing <- test(rnorm(10), rnorm(10))
    if(missing(parameter)) {
        ## Checking if the output is valid without a parameter value
        if(class(testing) == "numeric" && length(testing) == 1) {
            test_parameter <- TRUE
            param_element <- 1
        } else {
            if(any(names(testing) == "statistic")) {
                test_parameter <- TRUE
                param_element <- which(names(testing) == "statistic")
            } else {
                if(test.parameter) {
                    stop("Impossible to find a default parameter (statistic or simple output) to test.\nSpecify the parameter name with the argument:\nrand.test(..., parameter = \"my_parameter\")")
                } else {
                    warning("Impossible to find a default parameter (statistic or simple output) to output.\nplotting the results through plot(...) won't work directly.")
                }
            }
        }
    } else {
        param_element <- which(names(testing) == parameter)
        if(length(param_element) == 1) {
            test_parameter <- TRUE
        } else {
            stop("Impossible to find the parameter name in the test.")
        }
    }

    if(test.parameter) {
        ## p-value
        check.method(alternative, c("two-sided", "greater", "lesser"), msg = "alternative")
        
        ## Set p-value function
        if(alternative == "two-sided") {
            get.p.value <- function(random, observed, replicates) {
                ## Centring the randoms and observed
                center_random <- abs(random - mean(random))
                center_observed <- abs(mean(observed) - mean(random))
                ## Getting the p
                return((sum(center_random >= center_observed) + 1)/(replicates + 1))
            }
        }
        if(alternative == "greater") {
            get.p.value <- function(random, observed, replicates) {
                # Getting the p
                return((sum(random >= mean(observed)) + 1)/(replicates + 1))
            }
        }
        if(alternative == "lesser") {
            get.p.value <- function(random, observed, replicates) {
                # Getting the p
                return((sum(random <= mean(observed)) + 1)/(replicates + 1))
            }
        }
    }

    
    ## Select the sub-population size
    pop_size <- ifelse(is_rarefied, rarefaction, length(subset))

    ## Observed draw
    if(!is_rarefied) {
        observed_draws <- subset
        random_draws <- replicate(replicates, sample(1:length(distribution), length(subset)), simplify = FALSE)
    } else {
        observed_draws <- replicate(replicates, sample(subset, rarefaction), simplify = FALSE)
        random_draws <- replicate(replicates, sample(1:length(distribution), rarefaction), simplify = FALSE)
    }

    ## Applying test function
    if(resample) {
        lapply.test <- function(draw, test, distribution, ...) {
            return(test(distribution[draw], distribution, ...))
        }
    } else {
        lapply.test <- function(draw, test, distribution, ...) {
            return(test(distribution[draw], distribution[-draw], ...))
        }
    }


    ## Comparing the observed draws
    if(!is_rarefied) {
        observed_differences <- test(distribution[observed_draws], distribution)
    } else {
        observed_differences <- lapply(observed_draws, lapply.test, test, distribution, ...)
    }

    ## Comparing the random draws
    random_differences <- lapply(random_draws, lapply.test, test, distribution, ...)

    ## Unlist the random and observed differences
    random_parameters <- unlist(lapply(random_differences, function(X, param_element) X[[param_element]], param_element))
    if(!is_rarefied) {
        observed_parameters <- observed_differences[[param_element]]
    } else {
        observed_parameters <- unlist(lapply(observed_differences, function(X, param_element) X[[param_element]], param_element))
    }

    if(test.parameter) {
        ## Getting the test results
        if(!is_rarefied) {
            test_results <- c("Normal residuals" = (observed_parameters - mean(random_parameters)) / stats::sd(random_parameters), "Random mean" = mean(random_parameters), "Random variance" = stats::var(random_parameters))
        } else {
            test_results <- c("Mean Normal residuals" = mean((observed_parameters - mean(random_parameters)) / stats::sd(random_parameters)), "Random mean" = mean(random_parameters), "Random variance" = stats::var(random_parameters))
        }

        ## Calculating the p-value
        p_value <- get.p.value(random_parameters, observed_parameters, replicates)
    } 

    ## Making the results into a randtest object
    res <- list()

    ## Adding the default arguments
    res$rep <- replicates
    res$observed <- observed_differences
    res$random <- random_differences
    res$call <- match_call

    ## Adding the data
    res$sim <- random_parameters
    res$obs <- mean(observed_parameters)

    ## Adding the plot options (modified from ade4::as.randtest)
    r0 <- c(random_parameters, observed_parameters)
    l0 <- max(random_parameters, observed_parameters) - min(random_parameters, observed_parameters)
    w0 <- l0/(log(length(random_parameters), base = 2) + 1)
    xlim0 <- range(r0) + c(-w0, w0)
    h0 <- graphics::hist(random_parameters, plot = FALSE, nclass = 10)
    res$plot <- list(hist = h0, xlim = xlim0)

    ## Adding the test.parameter arguments
    res$alter <- ifelse(test.parameter, alternative, NA)
    res$pvalue <- ifelse(test.parameter, p_value, NA)
    if(test.parameter) {
        res$expvar <- test_results
    } else {
        res$expvar <- "No test was conducted on the results. Use plot(...) to visualise."
    }
        
    class(res) <- "randtest"

    return(res)
}