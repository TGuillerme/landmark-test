
# @param species the species name
# @param dataset the type of dataset (cranium or mandible)
# @param path the path where the processed data is

## Test pipeline
pipeline.test <- function(species, dataset, path, verbose = FALSE){
    ## Loading a dataset
    if(verbose) message("Load data...")
    load(paste0(path, species, ".Rda"))
    if(verbose) message("Done.\n")

    ## Selecting a dataset
    if(verbose) message("Select dataset...")
    data <- land_data[[dataset]]
    if(verbose) message("Done.\n")

    ## Procrustes variation ranges
    if(verbose) message("Calculate range...")
    procrustes_var <- variation.range(data$procrustes)
    if(verbose) message("Done.\n")

    ## landmarks partitions
    if(verbose) message("Determine partitions...")
    partitions <- list()
    for(part in unique(data$landmarkgroups[,2])) {
        partitions[[part]] <- which(data$landmarkgroups[,2] == part)
    }
    if(verbose) message("Done.\n")

    ## Size differences
    if(verbose) message("Run size difference test...")
    differences <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = area.diff)
    if(verbose) message("Done.\n")

    ## Probabilities of overlap
    if(verbose) message("Run overlap probability test...")
    overlaps <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = bhatt.coeff)
    if(verbose) message("Done.\n")

    return(list("differences" = differences, "overlaps" = overlaps, "species" = species, "dataset" = dataset))
}


## Summary pipeline
pipeline.plots <- function(results, export = FALSE){

    ## Plot the size differences
    make.plots(results$differences, type = "area difference", add.p = TRUE, correction = "bonferroni")

    ## Plot the size differences
    make.plots(results$overlaps, type = "Bhattacharyya Coefficient", add.p = TRUE, correction = "bonferroni")

}



## Function for applying the rand tests
lapply.rand.test <- function(partition, data, test, ...) {
    rand.test(data[, "radius"], partition, test = test, test.parameter = TRUE, ...)
}

## Function for plotting the test results
make.table <- function(results, correction) {

    ## Extract the values
    values <- lapply(results, function(X) return(c(X$obs, X$expvar, X$pvalue)))

    ## make the table
    summary_table <- cbind(seq(1:length(values)), do.call(rbind, values))

    ## Correcting the p.value
    if(!missing(correction)) {
        summary_table[, ncol(summary_table)] <- p.adjust(summary_table[, ncol(summary_table)],
            method = correction)
        p_name <- "p value (adjusted)"
    } else {
        p_name <- "p value"
    }

    ## Adding the colnames
    colnames(summary_table)[c(1,2,6)] <- c("Partition", "Observed", p_name)

    return(summary_table)
}


## Function for plotting the test results
make.plots <- function(results, type, add.p = FALSE, correction) {

    ## Number of plots
    n_plots <- length(results)

    ## Plotting parameters
    par(mfrow = c(ceiling(sqrt(n_plots)), floor(sqrt(n_plots))), bty = "n")

    ## Getting the results table
    if(add.p) {
        table_res <- make.table(results, correction)
    }

    for(one_plot in 1:n_plots) {
        ## Plot
        plot(results[[one_plot]], xlab = type, main = paste("Partition", one_plot))
        ## p_value
        if(add.p) {
            ## Get the coordinates for the text
            text_pos <- ifelse(table_res[one_plot, "Observed"] < table_res[one_plot, "Random mean"], "topleft", "topright")
            ## Add the text
            legend(text_pos, paste(colnames(table_res)[6], round(table_res[one_plot, 6], 5), sep = "\n")  , bty = "n")
        }
    }
}
