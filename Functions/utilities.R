# @param species the species name
# @param dataset the type of dataset (cranium or mandible)
# @param path the path where the processed data is
# @param combine.land whether to combine the landmarks (only two partitions)

## Test pipeline
pipeline.test <- function(species, dataset, path, verbose = FALSE, rarefaction, combine.land = FALSE, CI){
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
    procrustes_var <- variation.range(data$procrustes, type = "spherical", what = "radius", CI = CI)
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
        differences <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = area.diff, replicates = 1000, rarefaction = part_min, resample = FALSE)
    } else {
        differences <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = area.diff, replicates = 1000, resample = FALSE)
    }
    if(verbose) message("Done.\n")

    ## Probabilities of overlap
    if(verbose) message("Run overlap probability test...")
    if(rarefaction) {
        overlaps <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = bhatt.coeff, replicates = 1000, rarefaction = part_min, resample = FALSE)
    } else {
        overlaps <- lapply(partitions, lapply.rand.test, data = procrustes_var, test = bhatt.coeff, replicates = 1000, resample = FALSE)
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

    ## Rounding
    results_table[, c(2,3,6)] <- round(results_table[, c(2,3,6)], digits = digits)

    if(all(as.vector(round(results_table[, c(4,5)], digits = digits) == 0))) {
        results_table[, c(4,5)] <- round(results_table[, c(4,5)], digits = digits+2)
    } else {
        results_table[, c(4,5)] <- round(results_table[, c(4,5)], digits = digits)
    }

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
    textable <- xtable(results_table, caption = caption, label = label)

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

make.plots <- function(results, type, add.p = FALSE, correction, rarefaction = FALSE, rare.level, path) {

    ## Number of plots
    n_plots <- length(results)

    ## Plotting parameters
    if(!missing(path)) {
        pdf(file = path)
    } else {

        n_rows <- c(ceiling(sqrt(n_plots)), floor(sqrt(n_plots)))
        if(n_rows[1]*n_rows[2] < n_plots) {
            n_rows <- c(ceiling(sqrt(n_plots)), ceiling(sqrt(n_plots)))
        }

        par(mfrow = n_rows, bty = "n")
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

        ## Rarefaction (unless the rarefied results are invariant, i.e. minimum level)
        if(rarefaction && length(unique(unlist(results[[one_plot]]$observed))) != 1) {
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
#@param PointSize is for changing the size of spheres plotted

#Defining partitions using the define.module (needs individual execution)
plot.partitions<-function(land_data_partition, PartNames, PointSize){
  ##the object with the landmarks subset according to partitions
  Part=list()
  WomCrGPA<-land_data_partition$procrustes
  WomCrPart<-land_data_partition$landmarkgroups
  WomCrRef <- mshape(land_data_partition$procrustes$coords) 
  
  
  #provides the numbers of the parts
  PartLevels= unique(WomCrPart[,2])
  Colours<-rainbow(length(PartLevels))
  Colours <- c("blue", "orange", "green")
  
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
    spheres3d(WomCrRef[Part[[i]],1], WomCrRef[Part[[i]],2], WomCrRef[Part[[i]],3], col=Colours[i], lit=TRUE,radius = Pointsize, asp=F)
    
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
}
  
  
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



#@param CI: the confidence interval level
#@param rarefaction: whether to use the rarefied results (TRUE) or not (FALSE)
#@param print.token: whether to add the significance tokens (i.e. stars)
#@param rounding: the number of digits to print after 0.
#@param path: the path to the results
#@param species: the list of species as written in the results
#@param datasets: the names of the datasets as written in the results
#@param partitions: optional, the name of the landmark partitions
summarise.results <- function(CI, rarefaction, print.token = FALSE, rounding = 3, path = "../Data/Results/", species = c("Wombat", "Wombat_krefftii", "Wombat_latifrons", "Wombat_ursinus", "Wombat_lasiorhinus"), datasets = c("cranium", "mandible"), partitions = c("zygom", "restC", "snout", "back", "front", "restM")) {

    ## Printing significance tokens
    get.token <- function(p) {
        if(p > 0.05) {
            return("")
        }
        if(p < 0.05 && p > 0.01) {
            return(".")
        }
        if(p < 0.01 && p > 0.001) {
            return("*")
        }
        if(p < 0.001) {
            return("**")
        }
    }

    ## Make the empty table
    results_table <- data.frame(matrix(NA, ncol = length(species)*4, nrow = 6+1))

    ## Loop through the datasets
    for(sp in 1:length(species)) {

        for(ds in 1:length(datasets)) {

            ## Extract the results
            if(rarefaction) {
                load(paste0(path, species[sp], "_", datasets[ds], "rarefied_CI", CI, ".Rda"))
            } else {
                load(paste0(path, species[sp], "_", datasets[ds], "_CI", CI, ".Rda"))
            }

            ## Summarise the results
            difference <- make.table(results$difference, correction = "bonferroni")
            overlap <- make.table(results$overlaps, correction = "bonferroni")

            ## Fill the table
            if(ds == 1) {
                ## Values
                results_table[2:4, 1+(4*(sp-1))] <- round(difference[,2], digits = rounding)
                results_table[2:4, 3+(4*(sp-1))] <- round(overlap[,2], digits = rounding-1)
                ## Signif
                if(print.token) {
                    results_table[2:4, 2+(4*(sp-1))] <- paste0(round(difference[,6], digits = rounding),  sapply(difference[,6], get.token))
                    results_table[2:4, 4+(4*(sp-1))] <- paste0(round(overlap[,6], digits = rounding),  sapply(overlap[,6], get.token))
                } else {
                    results_table[2:4, 2+(4*(sp-1))] <- paste0(round(difference[,6], digits = rounding))
                    results_table[2:4, 4+(4*(sp-1))] <- paste0(round(overlap[,6], digits = rounding))
                }

            } else {
                ## Values
                results_table[5:7, 1+(4*(sp-1))] <- round(difference[,2], digits = rounding)
                results_table[5:7, 3+(4*(sp-1))] <- round(overlap[,2], digits = rounding-1)
                ## Signif
                if(print.token) {
                    results_table[5:7, 2+(4*(sp-1))] <- paste0(round(difference[,6], digits = rounding),  sapply(difference[,6], get.token))
                    results_table[5:7, 4+(4*(sp-1))] <- paste0(round(overlap[,6], digits = rounding),  sapply(overlap[,6], get.token))
                } else {
                    results_table[5:7, 2+(4*(sp-1))] <- paste0(round(difference[,6], digits = rounding))
                    results_table[5:7, 4+(4*(sp-1))] <- paste0(round(overlap[,6], digits = rounding))
                }
            }
        }
    }

    ## Renaming the table elements
    if(!missing(partitions)) {
        rownames(results_table) <- c("test", partitions)
    } else {
        rownames(results_table) <- c("test", paste("Cranium", c(1,2,3)), paste("Mandible", c(1,2,3)))
    }
    results_table[1,] <- rep(c("diff", "p", "overlap", "p"), length(species))
    colnames(results_table) <- names(unlist(sapply(species, rep, 4, simplify = FALSE)))

    ## Flip the table
    return(t(results_table))

}



# plot.test.results <- function(data, rarefaction, cols = c("grey", "green", "magenta"), ...) {

#     ## Returning significance
#     get.token <- function(p) {
#         if(p > 0.05) {
#             return(FALSE)
#         }
#         if(p < 0.05 && p > 0.01) {
#             return(FALSE)
#         }
#         if(p < 0.01 && p > 0.001) {
#             return(TRUE)
#         }
#         if(p < 0.001) {
#             return(TRUE)
#         }
#     }

# CI <- 95

# species <- c("Wombat", "Wombat_krefftii", "Wombat_latifrons", "Wombat_ursinus", "Wombat_lasiorhinus")
# datasets <- c("cranium", "mandible")

# rounding = 3

# results_table <- data.frame(matrix(NA, ncol = length(species)*4, nrow = 6+1))

# for(sp in 1:length(species)) {

#     for(ds in 1:length(datasets)) {

#         ## Extract the results
#         load(paste0("../Data/Results/", species[sp], "_", datasets[ds], "_CI", CI, ".Rda"))

#         ## Summarise the results
#         difference <- make.table(results$difference, correction = "bonferroni")
#         overlap <- make.table(results$overlaps, correction = "bonferroni")

#         ## Fill the table
#         if(ds == 1) {
#             ## Values
#             results_table[2:4, 1+(4*(sp-1))] <- round(difference[,2], digits = rounding)
#             results_table[2:4, 3+(4*(sp-1))] <- round(overlap[,2], digits = rounding-1)
#             ## Signif
#             results_table[2:4, 2+(4*(sp-1))] <- paste0(round(difference[,6], digits = rounding),  sapply(difference[,6], get.token)) # Remove the first component of the paste if only tokens needed
#             results_table[2:4, 4+(4*(sp-1))] <- paste0(round(overlap[,6], digits = rounding),  sapply(overlap[,6], get.token))

#         } else {
#             ## Values
#             results_table[5:7, 1+(4*(sp-1))] <- round(difference[,2], digits = rounding)
#             results_table[5:7, 3+(4*(sp-1))] <- round(overlap[,2], digits = rounding-1)
#             ## Signif
#             results_table[5:7, 2+(4*(sp-1))] <- paste0(round(difference[,6], digits = rounding),  sapply(difference[,6], get.token))
#             results_table[5:7, 4+(4*(sp-1))] <- paste0(round(overlap[,6], digits = rounding),  sapply(overlap[,6], get.token))
#         }
#     }
# }

# ## Renaming the table elements
# rownames(results_table) <- c("test", paste("Cranium", c(1,2,3)), paste("Mandible", c(1,2,3)))
# results_table[1,] <- rep(c("diff", "p", "overlap", "p"), length(species))
# colnames(results_table) <- names(unlist(sapply(species, rep, 4, simplify = FALSE)))

# ## Flip the table
# results_table <- t(results_table)


# ##
# kable(results_table, digits = 4)




#     ## Formatting data data as a matrix for image
#     age_rows <- list(c(1:7), c(8:14), c(15:21))
#     epoch_rows <- list(c(22:28), c(29:35), c(36:42))
#     signif_columns <- c(10:12)

#     ## Combining the age matrix
#     if(type == "age") {
#         data_matrix <- as.matrix(do.call(cbind, lapply(age_rows, function(rows, cols, data) return(data[rows, cols]), cols = signif_columns, data = data)))
#     } else {
#         data_matrix <- as.matrix(do.call(cbind, lapply(epoch_rows, function(rows, cols, data) return(data[rows, cols]), cols = signif_columns, data = data)))
#     }

#     ## Convert into binary
#     data_matrix <- ifelse(data_matrix, 1, 0)
#     data_matrix <- ifelse(is.na(data_matrix), 0.5, data_matrix)

#     ## Colours equivalent
#     request_colors <- length(unique(as.vector(data_matrix)))
#     if(request_colors == 3) {
#         ## Reorder colours (damn image()!)
#         cols <- cols[c(2,3,1)]
#     }

#     colours_table <- matrix(c(1,0,0.5), ncol = 1, dimnames = list(cols))
#     num_cols <- na.omit(sort(match(unique(as.vector(data_matrix)), colours_table)))

#     ## Selecting the colours
#     colours <- rownames(colours_table)[num_cols]

#     ## Plot the matrix
#     image(t(data_matrix[7:1,]), col = rev(colours), xaxt = "n", yaxt = "n", main = main, ...)
#     # image(t(data_matrix[7:1,]), col = rev(cols), xaxt = "n", yaxt = "n", main = main)
#     ## Add the lines
#     abline(v = c(0.3125, 0.6875), lty = 2)
#     ## Add the y axis
#     if(yaxis) axis(2, at = seq(from = 0, to = 1, by = 1/6), las = 2, label = rev(data$model[age_rows[[1]]]), tick = FALSE)
#     ## Add the x axis
#     if(xaxis) axis(1, at = c(0.125, 0.5, 0.85), label = c("stratigraphy", "duration", "number"), tick = FALSE)
#     ## Add the upper x axis
#     if(xaxis2) axis(3, at = c(seq(from = 0, to = 1, by = 1/8)), label = rep(c("e:1", "e:2", "e:3"), 3), tick = FALSE, padj = 1.5, cex = 0.8)
#     ## Add the right y axis
#     if(yaxis2) axis(4, at = 0.5, tick = FALSE, labels = type)
# }


# multi.plot.extinction <- function(data, type, cols = c("blue", "orange", "white"), main = "", data.names, ...) {

#     ## Check how many types
#     if(length(type) == 1) {

#         ## Setting up plot layout
#         plot_layout <- layout(matrix(c(1:length(data)), 1, length(data), byrow = TRUE), rep(1, length(data)), rep(1, length(data)), FALSE)
#         #layout.show(plot_layout)

#         ## First plot
#         par(mar = c(4, 6, 4, 0)) #c(bottom, left, top, right)
#         plot.extinction(data[[1]], type = type, col = cols, xaxis = TRUE, yaxis = TRUE, xaxis2 = TRUE, main = data.names[[1]], ...)

#         ## Other plots
#         for(slug in 2:length(data)) {
#             par(mar = c(4, 0, 4, 0)) #c(bottom, left, top, right)
#             plot.extinction(data[[slug]], type = type, col = cols, main = data.names[[slug]], xaxis = TRUE, xaxis2 = TRUE, ...)
#         }
#     } else {

#         ## Setting up plot layout
#         plot_layout <- layout(matrix(c(1:(length(data)*2)), 2, length(data), byrow = TRUE), rep(1, length(data)*2), rep(1, length(data)*2), FALSE)
#         #layout.show(plot_layout)

#         ## First plot (first row)
#         par(mar = c(0, 6, 4, 0)) #c(bottom, left, top, right)
#         plot.extinction(data[[1]], type = type[1], col = cols, yaxis = TRUE, xaxis2 = TRUE, main = data.names[[1]], ...)

#         # ## Other plots (first row)
#         # for(slug in 2:(length(data)-1)) {
#         #     par(mar = c(0, 0, 4, 0)) #c(bottom, left, top, right)
#         #     plot.extinction(data[[slug]], type = type[1], col = cols, main = data.names[[slug]], xaxis2 = TRUE, ...)
#         # }

#         ## Last plot (first row)
#         par(mar = c(0, 0, 4, 3)) #c(bottom, left, top, right)
#         plot.extinction(data[[length(data)]], type = type[1], col = cols, main = data.names[[length(data)]], xaxis2 = TRUE, yaxis2 = TRUE, ...)

#         ## First plot (second row)
#         par(mar = c(4, 6, 0, 0)) #c(bottom, left, top, right)
#         plot.extinction(data[[1]], type = type[2], col = cols, yaxis = TRUE, xaxis = TRUE, ...)

#         # ## Other plots (second row)
#         # for(slug in 2:(length(data)-1)) {
#         #     par(mar = c(4, 0, 0, 0)) #c(bottom, left, top, right)
#         #     plot.extinction(data[[slug]], type = type[2], col = cols, xaxis = TRUE, ...)
#         # }

#         ## Last plot (first row)
#         par(mar = c(4, 0, 0, 3)) #c(bottom, left, top, right)
#         plot.extinction(data[[length(data)]], type = type[2], col = cols, xaxis = TRUE, yaxis2 = TRUE, ...)
#     }
# }








