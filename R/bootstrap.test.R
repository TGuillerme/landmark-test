#' @title Bootstrap test
#'
#' @description Performs a bootstrap test
#'
#' @param distribution A \code{numeric} distribution to draw from.
#' @param subset A \code{numeric} vector of the elements of the distribution to test.
#' @param statistic The statistic to measure from the distribution (a \code{function}).
#' @param replicates A \code{numeric} value for the number of replicates (\code{default = 100}).
#' @param rarefaction \code{logical}, whether to reduce the size of the distribution to be equal to the subset (\code{TRUE}) or not (\code{FALSE}; default).
#' @param alternative Optional, if the parameter is tested, what is the alternative hypothesis. Can be \code{"two-sided"} (default), \code{"greater"} or \code{"lesser"}.
#' @param abs \code{logical}, whether to use absolute difference (\code{TRUE}) or not (\code{FALSE}; default).
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
#' ## Testing whether this partition has a different median than expected by change
#' boot_test <- bootstrap.test(var_range[, "radius"], random_part, statistic = median)
#'
#' ## Summarising the results
#' boot_test
#' 
#' ## Plotting the results
#' plot(boot_test)
#' 
#' @seealso \code{link{rand.test}}, \code{\link[ade4]{randtest}}
#' 
#' @author Thomas Guillerme
#' 
#' @export
#' @importFrom stats sd var
#' @importFrom graphics hist

bootstrap.test <- function(distribution, subset, statistic = mean, replicates = 100, rarefaction = FALSE, alternative = "two-sided", abs = FALSE, ...) {


    match_call <- match.call()

    ## Sanitizing
    ## Distribution and subset
    check.class(distribution, "numeric")
    check.class(subset, c("numeric", "integer"))
    if(length(subset) > length(distribution)) {
        stop("The subset is bigger than the whole distribution.")
    }

    ## Test
    check.class(statistic, "function")

    ## Replicates
    check.class(replicates, "numeric")
    check.length(replicates, 1, msg = " must be a single numeric value.")
    if(replicates < 1) {
        stop("At least one replicate must be run.")
    }

    ## p-value
    check.method(alternative, c("two-sided", "greater", "lesser"), msg = "alternative")

    ## abs
    check.class(abs, "logical")

    ## rarefaction
    check.class(rarefaction, "logical")    
    
    ## Set p-value function
    if(alternative == "two-sided") {
        get.p.value <- function(bs_t_statistic, t_statistic, replicates) {
            ## Centring the randoms and observed
            center_bs_t_statistic <- abs(bs_t_statistic - mean(bs_t_statistic))
            center_t_statistic <- abs(mean(t_statistic) - mean(bs_t_statistic))
            ## Getting the p
            return((sum(center_bs_t_statistic >= center_t_statistic))/(replicates))
        }
    }

    if(alternative == "greater") {
        get.p.value <- function(bs_t_statistic, t_statistic, replicates) {
            # Getting the p
            return((sum(bs_t_statistic >= mean(t_statistic)))/(replicates))
        }
    }
    if(alternative == "lesser") {
        get.p.value <- function(bs_t_statistic, t_statistic, replicates) {
            # Getting the p
            return((sum(bs_t_statistic <= mean(t_statistic)))/(replicates))
        }
    }

    ## Calculate the observed statistic
    statistic_subset <- statistic(distribution[subset])
    if(!rarefaction) {
        statistic_subset <- statistic(distribution[subset])
        statistic_distri <- statistic(distribution[-subset])
        t_statistic <- ifelse(abs, abs(statistic_subset - statistic_distri), statistic_subset - statistic_distri)
    } else {
        statistic_subset <- statistic(distribution[subset])
        statistic_distri <- replicate(replicates, statistic(sample(distribution[-subset], length(subset))))
        if(abs){
            t_statistic <- abs(statistic_subset - statistic_distri)
        } else {
            t_statistic <- statistic_subset - statistic_distri
        }

    }

    ## Bootstrap the data
    boostrap_subset <- replicate(replicates, distribution[sample(1:length(distribution), length(subset))], simplify = FALSE)
    if(!rarefaction) {
        boostrap_distri <- replicate(replicates, distribution[sample(1:length(distribution), length(distribution)-length(subset))], simplify = FALSE)
    } else {
        boostrap_distri <- replicate(replicates, distribution[sample(1:length(distribution), length(subset))], simplify = FALSE)
    }

    ## Calculate the statistics
    bs_stat_subset <- unlist(lapply(boostrap_subset, statistic))
    bs_stat_distri <- unlist(lapply(boostrap_distri, statistic))
    if(abs) {
        bs_t_statistic <- abs(bs_stat_subset - bs_stat_distri)
    } else {
        bs_t_statistic <- bs_stat_subset - bs_stat_distri
    }

    ## Calculate the p-value
    p_value <- get.p.value(bs_t_statistic, t_statistic, replicates)

    ## Get the test results
    if(!rarefaction) {
        test_results <- c("Residuals" = (t_statistic - mean(bs_t_statistic)) / stats::sd(bs_t_statistic), "Bootstrap mean" = mean(bs_t_statistic), "Bootstrap variance" = stats::var(bs_t_statistic))
    } else {
        test_results <- c("Residuals" = mean(t_statistic - bs_t_statistic) / stats::sd(bs_t_statistic), "Bootstrap mean" = mean(bs_t_statistic), "Bootstrap variance" = stats::var(bs_t_statistic))
    }

    ## Making the results into a randtest object
    res <- list()

    ## Adding the default arguments
    res$rep <- replicates
    res$call <- match_call

    ## Adding the data
    res$sim <- bs_t_statistic
    if(!rarefaction) {
        res$obs <- t_statistic
    } else {
        res$obs <- mean(t_statistic)
        res$observed <- t_statistic
    }

    ## Adding the plot options (modified from ade4::as.randtest)
    r0 <- c(bs_t_statistic, t_statistic)
    l0 <- max(bs_t_statistic) - min(bs_t_statistic)
    w0 <- l0/(log(length(bs_t_statistic), base = 2) + 1)
    xlim0 <- range(r0) + c(-w0, w0)
    h0 <- graphics::hist(bs_t_statistic, plot = FALSE, nclass = 10)
    res$plot <- list(hist = h0, xlim = xlim0)

    ## Adding the test.parameter arguments
    res$alter <- alternative
    res$pvalue <- p_value
    res$expvar <- test_results

    class(res) <- "randtest"

    return(res)
}
