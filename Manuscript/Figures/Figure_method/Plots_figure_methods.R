## Loading the libraries (and installing if necessary)
if(!require(devtools)) install.packages("devtools")
if(!require(knitr)) install.packages("knitr")
if(!require(geomorph)) install.packages("geomorph")
if(!require(dispRity)) install.packages("dispRity");library(dispRity)
if(!require(landmarktest)) install_github("TGuillerme/landmark-test");library(landmarktest)
source("../../../Functions/utilities.R")
set.seed(42)


## Loading a dataset
load("../../../Data/Processed/wombat_ursinus.Rda")

## Selecting a partition
cranium <- land_data$cranium

## Procrustes variation ranges
procrustes_var <- variation.range(cranium$procrustes, return.ID = TRUE)

## Partitions
partitions <- list()
for(part in unique(cranium$landmarkgroups[,2])) {
    partitions[[part]] <- which(cranium$landmarkgroups[,2] == part)
}

## Getting the zygomatic partition
random_part <- partitions[[1]]

## Measuring the difference for the zygomatic area
random_area <- rand.test(procrustes_var$range[, "radius"], random_part, test = area.diff,
                         test.parameter = TRUE)

## Measuring the difference for a zygomatic overlap
random_overlap <- rand.test(procrustes_var$range[, "radius"], random_part, test = bhatt.coeff,
                            test.parameter = TRUE)

## Plotting the area difference
plot.area <- function(distribution, X, reference, col) {
    par(bty = "n")
    plot(NULL, ylim = c(0, 0.01), xlim = c(1, length(X)), xlab = "", ylab = "")

    # par(bty = "n")
    # plot(NULL, ylim = c(0, 0.01), xlim = c(1, length(observed_draws)), xaxt = 'n', yaxt = 'n', xlab = "n", ylab = "n")
    polygon(
        x = c(1:length(X), rev(1:length(X))),
        y = c(sort(distribution[X], decreasing = TRUE), sort(distribution[reference], decreasing = FALSE)),
        col = col[3], border = NA
        )

    ## Add the lines
    lines(sort(distribution[X], decreasing = TRUE),         col = col[1], lwd = 2)
    lines(sort(distribution[reference], decreasing = TRUE), col = col[2], lwd = 2)
}

## Plotting the overlap function
plot.overlap  <- function(distribution, X, reference, col){
    ## Get the distributions
    distri_1 <- density(distribution[X])
    distri_2 <- density(distribution[reference])

    ## Get the x coordinates
    global_x <- c(distri_1$x, distri_2$x)
    x = seq(from = min(global_x), to = max(global_x), length.out = length(distri_1$y))

    ## Get the y coordinates
    y1 = distri_1$y
    y2 = distri_2$y

    ## Empty plot
    par(bty = "n")
    plot(NULL, xlim = c(-0.001, 0.01), ylim = c(0, 350), xlab = "", ylab = "")

    ## Plot the polygon
    polygon(x, pmin(y1, y2), col = col[3], border = NA)

    ## Add the lines
    lines(x, y1, type = "l", lwd = 2, col = col[1])
    lines(x, y2, type = "l", lwd = 2, col = col[2])
}





## Plot the results
pdf(file = "results_area.pdf")
par(bty = "n")
plot(random_area, main = "results area")
dev.off()
pdf(file = "results_overlap.pdf")
par(bty = "n")
plot(random_overlap, main = "results overlap")
dev.off()


## Setting the parameters for plotting the test data

## The data
distribution <- procrustes_var$range[,1]

## The different draws (first is the reference)
observed_draws <- partitions[[1]]
random_draws <- replicate(4, sample(1:length(distribution), length(observed_draws)), simplify = FALSE)

## The colours parameters
colours <- c("purple", "green", "grey") #obs, reference, difference


## Observed difference
## Overlap
pdf(file = "overlap_obs.pdf")
plot.overlap(distribution, observed_draws, reference = random_draws[[1]], col = colours)
dev.off()
## Area
pdf(file = "area_obs.pdf")
plot.area(distribution, observed_draws, reference = random_draws[[1]], col = colours)
dev.off()

## Simulated differences
for(one_plot in 2:4) {

    ## Overlap
    pdf(file = paste0("overlap_simul_0", one_plot-1, ".pdf"))
    plot.overlap(distribution, random_draws[[one_plot]], reference = random_draws[[1]], col = colours)
    dev.off()

    ## Area
    pdf(file = paste0("area_simul_0", one_plot-1, ".pdf"))
    plot.area(distribution, random_draws[[one_plot]], reference = random_draws[[1]], col = colours)
    dev.off()
}