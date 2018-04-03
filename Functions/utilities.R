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
    procrustes_var <- variation.range(data$procrustes, type = "spherical", what = "radius")
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
    rand.test(data[, "radius"], partition, test = test, test.parameter = TRUE, ...)
}
lapply.bootstrap.test <- function(partition, data, statistic, ...) {
    bootstrap.test(data[, "radius"], partition, statistic = statistic, ...)
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

make.xtable <- function(results, correction, digits = 3, caption, label, longtable = FALSE, path) {

    ## Make the results table
    results_table <- make.table(results, correction = correction)

    ## Add significance values
    p_col <- grep("p value", colnames(results_table))
    if(length(p_col) > 0) {
        for(row in 1:nrow(results_table)) {
            if(as.numeric(results_table[row, p_col]) < 0.05) {
                results_table[row, p_col] <- paste0("BOLD", results_table[row, p_col])
            }
        }
    }

    ## Bold cells function
    bold.cells <- function(x) gsub('BOLD(.*)', paste0('\\\\textbf{\\1', '}'), x)

    ## convert into xtable format
    textable <- xtable(results_table, digit = digits, caption = caption, label = label)

    if(!missing(path)) {
        if(longtable == TRUE) {
            cat(print(textable, tabular.environment = 'longtable', floating = FALSE, include.rownames = FALSE, sanitize.text.function = bold.cells), file = paste0(path, label, ".tex"))
        } else {
            cat(print(textable, include.rownames = FALSE, sanitize.text.function = bold.cells), file = paste0(path, label, ".tex"))
        }
    }

    if(longtable == TRUE) {
        print(textable, tabular.environment = 'longtable', floating = FALSE, include.rownames = FALSE, sanitize.text.function = bold.cells)
    } else {
        print(textable, include.rownames = FALSE, sanitize.text.function = bold.cells)
    }
}


## Function for plotting the test results
make.plots <- function(results, type, add.p = FALSE, correction, rarefaction = FALSE, rare.level, path) {

    ## Number of plots
    n_plots <- length(results)

    ## Plotting parameters
    if(!missing(path)) {
        pdf(file = path)
    } else {
        par(mfrow = c(ceiling(sqrt(n_plots)), floor(sqrt(n_plots))), bty = "n")
    }

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
            text_pos <- ifelse(table_res[one_plot, "Observed"] < table_res[one_plot, "Random mean"], "topleft", "topright")
            ## Add the text
            legend(text_pos, paste(colnames(table_res)[6], round(table_res[one_plot, 6], 5), sep = "\n")  , bty = "n")
        }
    }

    if(!missing(path)) {
        dev.off()
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

#procD code (for procD.allometry and procD.lm) analysis

#@param formula: a formula object (e.g. coords ~ Csize)
#@param procrustes: the procrustes object (e.g. land_data$cranium$procrustes)
#@param procD.fun: the procD function (e.g. procD.allometry)
#@param ...: any optional arguments to be passed to procD.fun (e.g. logsz = FALSE, iter = 1, etc...)
handle.procD.formula <- function(formula, procrustes, procD.fun = procD.allometry, ...) {
  
  geomorph_data_frame <- geomorph.data.frame(procrustes)
  
  return(procD.fun(f1 = formula, data = geomorph_data_frame, ...))
}


#Allometry analysis (based on above handle.procD.formula)

allom.shape<-function (procrustes_coordinate_file_with_centroid_size){
  
  Allometry <- handle.procD.formula(formula=coords~ Csize, procrustes=procrustes_coordinate_file_with_centroid_size, procD.fun = procD.allometry, logsz = FALSE, iter = 1000)
  print(attributes(Allometry))
  return(Allometry)
  
  
  
#Reducing datasets to those with counterparts
#@params AllData is a list of procrustes objects after gpa (e.g. land_data, in this case the different species); AllClassifiers is a list of classifiers matched with the AllData shape dataset, which includes subsetting information  
  
reduce.check<-function(AllData, AllClassifiers){
    
    coords_PLS_output=list()
    check_output=list()
    
    for (i in 1:length(AllData)){
      coords_PLS_output[[i]]=list()
      
      for (k in 1:length(AllClassifiers[[i]])){
        coords_PLS_output[[i]][[k]] <- AllData [[i]][[k]]$procrustes$coords [ , ,as.character(AllClassifiers[[i]][[k]]$TBPLS) != "Nil"]
      }
      names(coords_PLS_output[[i]])<-names(AllData[[i]])
    }
    
    names(coords_PLS_output)<-names(AllData)
    
    for (i in 1:length(AllData)){
      matchCheck=match(attributes(coords_PLS_output[[i]]$cranium)$dimnames[[3]], attributes(coords_PLS_output[[i]]$mandible)$dimnames[[3]])
      check_output[[i]]=!is.na(matchCheck)&&all(matchCheck==sort(matchCheck))
    }
    names(check_output)<-names(AllData)
    return(list(coords_PLS_output, check_output))
  }
}