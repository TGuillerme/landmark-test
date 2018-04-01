# @param species the species name
# @param dataset the type of dataset (cranium or mandible)
# @param path the path where the processed data is
# @param combine.land whether to combine the landmarks (only two partitions)

## Test pipeline
pipeline.test <- function(species, dataset, path, verbose = FALSE, rarefaction, combine.land = FALSE){
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
    procrustes_var <- variation.range(data$procrustes, type = "vector", what = "length")
    if(verbose) message("Done.\n")

    ## landmarks partitions
    if(verbose) message("Determine partitions...")
    partitions <- list()

    ## Combine the landmarks
    if(combine.land) {
        ## Select the biggest partitions
        biggest_partition <- which(table(data$landmarkgroups[,2]) == max(table(data$landmarkgroups[,2])))
        ## Combine landmarks not in that partition
        data$landmarkgroups[,2] <- ifelse(data$landmarkgroups[,2] != biggest_partition, 1, biggest_partition)
    }

    for(part in 1:length(unique(data$landmarkgroups[,2]))) {
        partitions[[part]] <- which(data$landmarkgroups[,2] == unique(data$landmarkgroups[,2])[part])
    }

    if(missing(rarefaction)) {
        rarefaction <- FALSE
    }

    if(rarefaction) {
        part_min <- min(unlist(lapply(partitions, length)))
    }

    if(verbose) message("Done.\n")

    ## Size differences
    if(verbose) message("Run size difference test...")
    if(rarefaction) {
        differences <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = area.diff, replicates = 1000, rarefaction = part_min)
    } else {
        differences <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = area.diff, replicates = 1000)
    }
    if(verbose) message("Done.\n")

    ## Probabilities of overlap
    if(verbose) message("Run overlap probability test...")
    if(rarefaction) {
        overlaps <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = bhatt.coeff, replicates = 1000, rarefaction = part_min)
    } else {
        overlaps <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = bhatt.coeff, replicates = 1000)
    }
    if(verbose) message("Done.\n")

    if(!rarefaction) {
        return(list("differences" = differences, "overlaps" = overlaps, "species" = species, "dataset" = dataset))
    } else {
        return(list("differences" = differences, "overlaps" = overlaps, "species" = species, "dataset" = dataset, "rarefaction" = part_min))
    }
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
    rand.test(data[, "length"], partition, test = test, test.parameter = TRUE, ...)
}
lapply.bootstrap.test <- function(partition, data, statistic, ...) {
    bootstrap.test(data[, "length"], partition, statistic = statistic, ...)
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
make.plots <- function(results, type, add.p = FALSE, correction, rarefaction = FALSE, rare.level) {

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

        if(!rarefaction) {
            main_lab <- paste("Partition", one_plot)
        } else {
            main_lab <- paste0("Partition ", one_plot, " (rarefied - ", rare.level, ")")
        }

        plot(results[[one_plot]], xlab = type, main = main_lab)

        ## Rarefaction
        if(rarefaction) {
            add.rare.plot(results[[one_plot]])
        }

        ## p_value
        if(add.p) {
            ## Get the coordinates for the text
            text_pos <- ifelse(table_res[one_plot, "Observed"] < table_res[one_plot, 2], "topleft", "topright")
            ## Add the text
            legend(text_pos, paste(colnames(table_res)[6], round(table_res[one_plot, 6], 5), sep = "\n")  , bty = "n")
        }
    }
}


# Utilities based on existing functions

#colouring partition spheres
#@param land_data_partition the landmark data e.g. land_data$cranium
#@param partnames is an optional vector with names for each partition number

plot.partitions<-function(land_data_partition, PartNames){
  ##the object with the landmarks subset according to partitions
  Part=list()
  WomCrGPA<-land_data_partition$procrustes
  WomCrPart<-land_data_partition$landmarkgroups
  WomCrRef <- mshape(WomCrGPA$coords) 
  
  
  #provides the numbers of the parts
  PartLevels= unique(WomCrPart[,2])
  Colours<-rainbow(length(PartLevels))
  
  ##subset the landmarks according to the partitions
  for(i in 1:length(PartLevels)){
    Part[[i]]<-which (WomCrPart[,2] == PartLevels[[i]])
  }
  
  ##provides names for each of the partitions (optional and requires a name vector to be given)
  if (!missing(PartNames)){
    for (i in 1:length(PartLevels)){
      names(Part)[i]<-PartNames[i]
    }
  }
  ##colours the spheres for each partition
  open3d()
  for (i in 1:length(PartLevels)){
    spheres3d(WomCrRef[Part[[i]],1], WomCrRef[Part[[i]],2], WomCrRef[Part[[i]],3], col=Colours[i], lit=TRUE,radius = 0.001, asp=F)
    
  }
  
}

#Visualizing differences between PC min and max
#@param x is the coordinates after gpa e.g. land_data$cranium
#@param minfirst is whether min vs max or other way round (for PlotRefToTarget )

PCA.vectors<-function(x, minfirst=TRUE){
  gridPar=gridPar(pt.bg = "white", pt.size = 0.5)#These are needed to give it the right parameters for the point size, colour and size - the default is too large points
  open3d()
  
  PCA=plotTangentSpace(x$procrustes$coords)
  open3d()
  if (minfirst==TRUE){
    plotRefToTarget(PCA$pc.shapes$PC1min,PCA$pc.shapes$PC1max, method="vector", gridPars=gridPar, label = F)
  } else {
    plotRefToTarget(PCA$pc.shapes$PC1max,PCA$pc.shapes$PC1min, method="vector", gridPars=gridPar, label = F)
  }
}

#allometry analysis

#@param procrustes_coordinate_file_with_centroid_size: "land_data_cranium_procrustes" of land_data contains all the objects of gpagen, including coords and Csize

#this does not work because the code fails to extract the coords within the loop. but the workflow works outside the loop
allom.shape<-function (procrustes_coordinate_file_with_centroid_size){
  coords=procrustes_coordinate_file_with_centroid_size$coords
  Csize=procrustes_coordinate_file_with_centroid_size$Csize
  
  Allometry <- procD.allometry(coords~ Csize, f2=NULL, f3=NULL, logsz = FALSE, iter=999)
  return(Allometry)
}


