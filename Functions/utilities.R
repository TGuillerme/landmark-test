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


## Function for translating the names in the datasets into the actual species names
translate.name <- function(name) {
    names <- c("All species", "Lasiorhinus", "Lasiorhinus krefftii", "Lasiorhinus latifrons", "Vombatus ursinus")

    if(name == "Wombat_ursinus") {
        return(names[5])
    }
    if(name == "Wombat") {
        return(names[1])
    }
    if(name == "Wombat_lasiorhinus") {
        return(names[2])
    }
    if(name == "Wombat_krefftii") {
        return(names[3])
    }
    if(name == "Wombat_latifrons") {
        return(names[4])
    }
    return(name)
}

## Function for getting the partition order and names
get.partition.args <- function(partition) {
    if(partition == "cranium") {
        return(list("partitions" = c("Zygomatic Arch", "Tip of snout", "Remainder"), "partitions.order" = c(1, 3, 2)))
    }
    if(partition == "mandible") {
        return(list("partitions" = c("Masticatory insertions", "Symphyseal area", "Remainder"), "partitions.order" = c(3, 1, 2)))
    }
    return(list(NULL))
}

## Function for plotting the test results
make.table <- function(results, correction, partition.args) {

    if(missing(partition.args)) {
        partition.args <- list(NULL)
    }

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

    ## Adding the partitions names
    if(!is.null(partition.args$partitions)) {
        summary_table <- data.frame(summary_table)
        summary_table[,1] <- partition.args$partitions
    }

    ## Adding the colnames
    colnames(summary_table)[c(1,2,6)] <- c("Partition", "Observed", p_name)

    ## Re-ordering partitions (if not missing)
    if(!is.null(partition.args$order.partitions)) {
        summary_table <- summary_table[, partition.args$order.partitions]
    }

    return(summary_table)
}

make.xtable <- function(results, correction, partition.args, digits = 3, caption, label, longtable = FALSE, path) {

    ## Make the results table
    results_table <- make.table(results, correction = correction, partition.args)

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
    spheres3d(WomCrRef[Part[[i]],1], WomCrRef[Part[[i]],2], WomCrRef[Part[[i]],3], col=Colours[i], lit=TRUE,radius = PointSize, asp=F)
    
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

heatplot.PCs(CW$cranium, minfirst=FALSE, PC_axis=1)

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
#@param partitions.order: optional, reordering the partition names
#@param partitions: optional, the name of the landmark partitions columns
#@param species.names: the names of species to display
summarise.results <- function(CI, rarefaction, print.token = FALSE, rounding = 3, path = "../Data/Results/", species = c("Wombat", "Wombat_lasiorhinus", "Wombat_krefftii", "Wombat_latifrons", "Wombat_ursinus"), datasets = c("cranium", "mandible"), partitions.order = c(1, 3, 2, 6, 4, 5), partitions, species.names) {

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

    ## Reordering the rows (future columns)
    results_table[2:7,] <- results_table[(partitions.order+1),]

    ## Renaming the table elements
    if(!missing(partitions)) {
        rownames(results_table) <- c("test", partitions)
    } else {
        rownames(results_table) <- c("test", paste("Cranium", c(1,2,3)), paste("Mandible", c(1,2,3)))
    }
    results_table[1,] <- rep(c("diff", "p", "overlap", "p"), length(species))
    if(!missing(species.names)) {
        colnames(results_table) <- names(unlist(sapply(species.names, rep, 4, simplify = FALSE)))
    } else {
        colnames(results_table) <- names(unlist(sapply(species, rep, 4, simplify = FALSE)))
    }

    ## Flip the table
    return(t(results_table))

}

#@param data: the non-rarefied summarised data
#@param rarefaction: the rarefied summarised data
#@param no.rar: optional, which columns to highlight as not rarefied (e.g. cranium 3 and mandible 1: no.rar = c(3,4))
#@param ignore.non.signif: whether to ignore the non-significant results for rarefaction highlights (TRUE) or not (FALSE).
#@pram partitions: the names of the partitions (c("cranium", "mandible"))
#@param cols: a vector of three colours for each pixel, the first one is non-significant results, the second one is when the difference is significant but not the overlap and the third one is when everything is significant
#@param threshold: the significance threshold (default = 0.01)
#@param ylabs: the labels for the y axis ("all_data", "Vombatus", etc). If missing the ones from data are used.
#@param xlabs: the labels for the x axis ("cranium1", etc). If missing the ones from data are used.
plot.test.results <- function(data, rarefaction, no.rar, ignore.non.signif = TRUE, partitions = c(expression(bold("Cranium")), expression(bold("Mandible"))), cols = c("grey", "magenta", "green"), ylabs, xlabs) {

    ## Making the x and y labels (if needed)
    if(missing(ylabs)) {
        ylabs <- gsub("1", "", rownames(data)[seq(from = 1, to = nrow(data), by = 4)])
    }
    if(missing(xlabs)) {
        xlabs <- colnames(data[, -1])
    }

    ## Converting the matrix in to numeric blocks
    make.blocks <- function(col) {
        from <- as.list(seq(from = 1, to = length(col), by = 4))
        to <- as.list(seq(from = 0, to = length(col), by = 4)[-1])
        return(mapply(function(from, to, col) return(col[from:to]), from, to, MoreArgs = list(col = col), SIMPLIFY = FALSE))
    }
    blocks <- apply(apply(data[,-1], 2, as.numeric), 2, make.blocks)

    ## Selecting the threshold level
    level.selector <- function(block, threshold = 0.01) {
        return(ifelse(block[2] > threshold, 1, ifelse(block[4] > threshold, 2, 3)))
    }
    ## Transform the list of blocks in an image matrix
    image_matrix <- matrix(unlist(lapply(blocks, lapply, level.selector)), ncol = ncol(data[, -1]), byrow = FALSE)

    ## Plot the main image
    par(mar = c(2, max(nchar(ylabs))/2, 4, 2)) #c(bottom, left, top, right)
    #image(t(image_matrix[nrow(image_matrix):1,]), col = cols, xaxt = "n", yaxt = "n", ...)
    image(t(image_matrix[nrow(image_matrix):1,]), col = cols, xaxt = "n", yaxt = "n")

    ## Adding the y labels
    axis(2, at = seq(from = 0, to = 1, length.out = nrow(image_matrix)), las = 2, label = rev(ylabs), tick = FALSE)

    ## Add the x labels
    if(length(grep("\\n", xlabs) > 0)) {
        padj <- 0.5
    } else {
        padj <- 1
    }

    axis(3, at = seq(from = 0, to = 1, length.out = ncol(image_matrix)), label = xlabs, tick = FALSE, padj = padj)
    axis(3, at = c(0.25, 0.75), label = partitions, tick = FALSE, padj = -2)

    ## Adding the values
    value_to_plot <- ifelse(image_matrix != 1, TRUE, FALSE)
    rownames(value_to_plot) <- rev(seq(from = 0, to = 1, length.out = nrow(value_to_plot)))
    colnames(value_to_plot) <- seq(from = 0, to = 1, length.out = ncol(value_to_plot))

    ## Getting the list of blocks coordinates from a named TRUE/FALSE matrix
    get.coords <- function(list_coords) {
        list_blocks_x <- as.numeric(t(apply(list_coords, 1, function(x) names(x))))
        list_blocks_y <- as.numeric(apply(list_coords, 2, function(x) names(x)))
        coords_x <- list_blocks_x[list_coords]
        coords_y <- list_blocks_y[list_coords]
        return(list(coords_x, coords_y))
    }

    values_coords <- get.coords(value_to_plot)
    values <- unlist(lapply(blocks, lapply, function(x) return(x[1])))[value_to_plot]
    text(x = values_coords[[1]], y = values_coords[[2]], labels = values)

    ## Adding the rarefaction (square the pixel if equal)
    blocks <- apply(apply(rarefaction[,-1], 2, as.numeric), 2, make.blocks)
    ## Transform the list of blocks in an image matrix
    image_rar <- matrix(unlist(lapply(blocks, lapply, level.selector)), ncol = ncol(rarefaction[, -1]), byrow = FALSE)
    ## Getting the rarefaction coordinates
    if(ignore.non.signif) {
        rar_coords <- image_rar == ifelse(image_matrix == 1, 0, image_matrix)
    } else {
        rar_coords <- image_rar == image_matrix
    }
    rownames(rar_coords) <- rev(seq(from = 0, to = 1, length.out = nrow(image_matrix)))
    colnames(rar_coords) <- seq(from = 0, to = 1, length.out = ncol(image_matrix))

    ## Removing no.rar columns if not missing
    if(!missing(no.rar)) {
        rar_coords[,no.rar] <- FALSE
    }

    ## Getting the list of blocks coordinates that are the same between rar and normal
    rar_coords <- get.coords(rar_coords)

    ## Getting the polygon coordinates
    get.polygon.coordinates <- function(block_x, block_y, image_matrix, lwd) {
        ## Get the size of the pixels
        x_size <- diff(seq(from = 0, to = 1, length.out = ncol(image_matrix)))[1]
        y_size <- diff(seq(from = 0, to = 1, length.out = nrow(image_matrix)))[1]

        if(!missing(lwd)) {
            x_size <- x_size - lwd/1000
            y_size <- y_size - lwd/1000
        }

        x_coords <- c(block_x-x_size/2, block_x-x_size/2, block_x+x_size/2, block_x+x_size/2)
        y_coords <- c(block_y-y_size/2, block_y+y_size/2, block_y+y_size/2, block_y-y_size/2)

        return(list(x_coords, y_coords))
    }

    ## Adding the polygons
    for(poly in 1:length(rar_coords[[1]])) {
        block_polygon <- get.polygon.coordinates(rar_coords[[1]][poly], rar_coords[[2]][poly], image_matrix)
        polygon(x = block_polygon[[1]], y = block_polygon[[2]], lwd = 3)
    }

    ## Separator
    abline(v = 0.5, lty = 2)
}

