#' @title Plot shape differences as geomorph::PlotRefToTarget with heat
#'
#' @description Modification of geomorph::PlotRefToTarget to allow to add heat maps (only works for \code{methods = c("points", "vector")})
#'
#' @param M1 Matrix of landmark coordinates for the first (reference) specimen
#' @param M2 Matrix of landmark coordinates for the second (target) specimen
# @param mesh A mesh3d object for use with {method="surface"}
# @param outline An x, y curve or curves warped to the reference (2D only)
#' @param method Method used to visualize shape difference; see below for details
# @param dots$mag The desired dots$magnification to be used when visualizing the shape difference (e.g., dots$mag=2)
# @param dots$links An optional matrix defining for dots$links between landmarks
# @param dots$label A logical value indicating whether landmark numbers will be plotted
# @param dots$axes A logical value indicating whether the box and dots$axes should be plotted (points and vector only)
#' @param gridPars An optional object made by \code{\link{gridPar}}
# @param useRefPts An option (logical value) to use reference configuration points rather than target configuration points (when {method = "TPS"})
#' @param ... Additional parameters not covered by \code{\link{gridPar}} to be passed to \code{\link{plot}}, \code{\link{rgl::plot3d}} or \code{\link{shade3d}}
#' @param col.points A single colour value, a vector of colour values or a colour \code{function} to colour the points (if \code{method = "points"})
#' @param val.points A \code{vector} of length \code{dim(M1)[1]} to attribute the points colour values (if \code{method = "points"})
#' @param col.vector A single colour value, a vector of colour values or a colour \code{function} to colour the vectors (if \code{method = "vector"})
#' @param val.vector A \code{vector} of length \code{dim(M1)[1]} to attribute the vectors colour values (if \code{method = "vector"})
#'
#' @export
#' @author Dean Adams, Emma Sherratt & Michael Collyer - modified by Thomas Guillerme
#' 
#' 

# library(geomorph)
# data(scallops)
# procrustes <- gpagen(scallops$coorddata)
# variation <- variation.range(procrustes, return.ID = TRUE)
# procrustes_var <- variation$range
# min_max <- variation$min.max

# M1 <- procrustes$coords[, , min_max[2]]
# M2 <- procrustes$coords[, , min_max[1]]

# method <- "vector"
# gridPars = gridPar(pt.bg = "white", pt.size = 0.5)

# col.points <- heat.colors
# val.points <- procrustes_var[, 1]
# col.vector <- heat.colors
# val.vector <- procrustes_var[, 1]


# gridPars = gridPar(pt.bg = "white", pt.size = 0.5)
# open3d()
# plotRefToTarget.heat(M1, M2, method = "vector", griPars = gridPar, col.points = heat.colors, val.points = val.points, col.vector = heat.colors, val.vector = val.vector)

plotRefToTarget.heat <- function(M1, M2, method = c("vector", "points"), gridPars = NULL, col.points, val.points, col.vector, val.vector, ...) {

  ## Get the methods (defaults)
  method <- match.arg(method)

  ## Manage default arguments from plotRefToTarget
  dots <- list(...)

  ## dots$mag
  if(is.null(dots$mag)) {
    dots$mag <- 1.0
  }

  ## dots$label
  if(is.null(dots$label)) {
    dots$label <- FALSE
  }

  ## dots$axes
  if(is.null(dots$axes)) {
    dots$axes <- FALSE
  }

  ## useRefPts
  if(is.null(dots$useRefPts)) {
    dots$useRefPts <- FALSE
  }

  ## Check M1 and M2
  if(any(is.na(M1)) == TRUE) {
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if(any(is.na(M2)) == TRUE) {
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")
  }
 
  ## Check griPars
  if(is.null(gridPars)) gP = gridPar() else gP=gridPars

  k <- dim(M1)[2]

  dots$mag <- (dots$mag - 1)
  M2 <- M2 + (M2 - M1) * dots$mag
  
  limits = function(x, s) { 
      r = range(x)
      rc = scale(r, scale = FALSE)
      l = mean(r) + s * rc
  }

  ## Colour gradient function
  get.colour.gradient <- function(values, col.fun) {

    ## Sort the data by range
    histo <- hist(values, plot = FALSE)
    n_col <- length(histo$counts)

    ## Get the gradient
    avail_cols <- rev(col.fun(n_col))

    ## Attribute the colours
    col_out <- rep(avail_cols[n_col], length(values))

    ## Attribute the colour by range
    for(colour in rev(1:n_col)) {
      col_out[(values <= histo$breaks[colour])] <- avail_cols[colour]
    }

    return(col_out)
  }

  ## Points colours
  if(missing(col.points)) {
    col_points <- "black"
  } else {
    if(class(col.points) != "function") {
      col_points <- col.points
    } else {
      if(!missing(val.points)) {
        col_points <- get.colour.gradient(val.points, col.points)
      } else {
        stop("val.points is missing to apply the col.points function.")
      }
    }
  }

  ## Vector colours
  if(missing(col.vector)) {
    col_vector <- "black"
  } else {
    if(class(col.vector) != "function") {
      col_vector <- col.vector
    } else {
      if(!missing(val.points)) {
      col_vector <- get.colour.gradient(val.vector, col.vector)
      } else {
        stop("val.points is missing to apply the col.points function.")
      }
    }
  }

  ## 2D plots
  if(k == 2) {

        if(method == "vector") {
            if(dots$axes == TRUE) {
              plot(M1, asp = 1, type = "n", xlab = "x", ylab = "y", xlim = limits(M1[, 1], 1.25), ylim = limits(M1[, 2], 1.25), ...)
            } else {
              plot(M1, asp = 1, type = "n", xlab = "", ylab = "", xlim = limits(M1[, 1], 1.25), axes = FALSE, ylim = limits(M1[, 2], 1.25), ...)
            }
            if(is.null(dots$links) == FALSE) {
                linkcol <- rep(gP$link.col, nrow(dots$links))[1:nrow(dots$links)]
                linklwd <- rep(gP$link.lwd, nrow(dots$links))[1:nrow(dots$links)]
                linklty <- rep(gP$link.lty, nrow(dots$links))[1:nrow(dots$links)]
                
                for (i in 1:nrow(dots$links)) {
                    segments(M2[dots$links[i, 1], 1], M2[dots$links[i, 1], 2], M2[dots$links[i, 2], 1], M2[dots$links[i, 2], 2], col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
                }
            }
            if(dots$label == TRUE) {
              text(M1, label = paste(1:dim(M1)[1]), adj = gP$txt.adj, pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)
            }

            arrows(M1[, 1], M1[, 2], M2[, 1], M2[, 2], length = 0.075, lwd = 2, col = col_vector)

            points(M1, pch = 21, bg = col_points, cex = gP$pt.size, col = gP$pt.bg)
        }

        if(method == "points") {
            if(dots$axes == TRUE) {
              plot(M1, asp = 1, pch = 21, type = "n", xlim = limits(M1[, 1], 1.25), ylim = limits(M1[, 2], 1.25), xlab = "x", ylab = "y", ...)
            }
            if(dots$axes == FALSE) {
                plot(M1, asp = 1, pch = 21, type = "n", xlim = limits(M1[, 1], 1.25), axes = FALSE, ylim = limits(M1[, 2], 1.25), xlab = "", ylab = "", ...)
            }
            if(dots$label == TRUE) {text(M1, label = paste(1:dim(M1)[1]), adj = gP$txt.adj, pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)}
            if(!is.null(outline)) {
                curve.warp <- geomorph::tps2d(outline, M1, M2)
                points(outline, pch = 19, cex = gP$out.cex, col = gP$out.col) 
                points(curve.warp, pch = 19, cex = gP$tar.out.cex, col = gP$tar.out.col) 
            }
            if(is.null(dots$links) == FALSE) {
                linkcol <- rep(gP$link.col, nrow(dots$links))[1:nrow(dots$links)]
                linklwd <- rep(gP$link.lwd, nrow(dots$links))[1:nrow(dots$links)]
                linklty <- rep(gP$link.lty, nrow(dots$links))[1:nrow(dots$links)]
                tarlinkcol <- rep(gP$tar.link.col, nrow(dots$links))[1:nrow(dots$links)]
                tarlinklwd <- rep(gP$tar.link.lwd, nrow(dots$links))[1:nrow(dots$links)]
                tarlinklty <- rep(gP$tar.link.lty, nrow(dots$links))[1:nrow(dots$links)]
                for (i in 1:nrow(dots$links)) {
                    segments(M1[dots$links[i, 1], 1], M1[dots$links[i, 1], 2], M1[dots$links[i, 2], 1], M1[dots$links[i, 2], 2], col = linkcol[i],lty = linklty[i], lwd = linklwd[i])
                    segments(M2[dots$links[i, 1], 1], M2[dots$links[i, 1], 2], M2[dots$links[i, 2], 1], M2[dots$links[i, 2], 2], col = tarlinkcol[i], lty = tarlinklty[i], lwd = tarlinklwd[i])
                }
            }
            points(M2, pch = 21, bg = col_points, cex = gP$tar.pt.size)
            points(M1, pch = 21, bg = col_points, cex = gP$pt.size)
        }         
    }

    ## 3D plots
    if(k == 3) {

        if(method == "vector") {

            if(dots$axes == TRUE) {
              rgl::plot3d(M1, type = "s", col = col_vector, size = gP$pt.size, aspect = FALSE, ...)
            }
            if(dots$axes == FALSE) {
                rgl::plot3d(M1, type = "s", col = col_vector, size = gP$pt.size, aspect = FALSE, xlab = "", ylab = "", zlab = "", axes = F, ...)
            }
            if(dots$label == TRUE) {
              rgl::text3d(M1, texts = paste(1:dim(M1)[1]), adj = (gP$txt.adj + gP$pt.size), pos = (gP$txt.pos + gP$pt.size), cex = gP$txt.cex, col = gP$txt.col)
            }
            
            

            ## Use a mapply here!!!!!

            add.segment <- function(i, M1, M2, lwd = 2, col) {
              rgl::segments3d(rbind(M1[i, ], M2[i, ]), lwd = lwd, col = col[i])
            }

            lapply(as.list(1:nrow(M1)), add.segment, M1 = M1, M2 = M2, col = col_vector)

            # for (i in 1:nrow(M1)) {
            #     rgl::segments3d(rbind(M1[i, ], M2[i, ]), lwd = 2, col = col_vector[i])
            # }




            if(is.null(dots$links) == FALSE) {
                tarlinkcol <- rep(gP$tar.link.col, nrow(dots$links))[1:nrow(dots$links)]
                tarlinklwd <- rep(gP$tar.link.lwd, nrow(dots$links))[1:nrow(dots$links)]
                tarlinklty <- rep(gP$tar.link.lty, nrow(dots$links))[1:nrow(dots$links)]
                for (i in 1:nrow(dots$links)) {
                    rgl::segments3d(rbind(M2[dots$links[i, 1], ], M2[dots$links[i, 2], ]),  col = tarlinkcol[i], lty = tarlinklty[i], lwd = tarlinklwd[i])
                }
            }
        }
        if(method == "points") {
            if(dots$axes == TRUE) {
              rgl::plot3d(M1, type = "s", col = col_points, size = gP$pt.size, aspect = FALSE, ...)
              rgl::plot3d(M2, type = "s", col = gP$tar.pt.bg, size = gP$tar.pt.size, add = TRUE)
            }
            if(dots$axes == FALSE) {
                rgl::plot3d(M1, type = "s", col = col_points, size = gP$pt.size, aspect = FALSE, xlab = "", ylab = "", zlab = "", axes = F, ...)
                rgl::plot3d(M2, type = "s", col = gP$tar.pt.bg, size = gP$tar.pt.size, add = TRUE)
            }
            
            if(dots$label == TRUE) {
              rgl::text3d(M1, texts = paste(1:dim(M1)[1]), adj = (gP$txt.adj+gP$pt.size), pos = (gP$txt.pos+gP$pt.size), cex = gP$txt.cex, col = gP$txt.col)
            }
            if(is.null(dots$links) == FALSE) {
                linkcol <- rep(gP$link.col, nrow(dots$links))[1:nrow(dots$links)]
                linklwd <- rep(gP$link.lwd, nrow(dots$links))[1:nrow(dots$links)]
                linklty <- rep(gP$link.lty, nrow(dots$links))[1:nrow(dots$links)]
                tarlinkcol <- rep(gP$tar.link.col, nrow(dots$links))[1:nrow(dots$links)]
                tarlinklwd <- rep(gP$tar.link.lwd, nrow(dots$links))[1:nrow(dots$links)]
                tarlinklty <- rep(gP$tar.link.lty, nrow(dots$links))[1:nrow(dots$links)]
                
                for (i in 1:nrow(dots$links)) {
                    rgl::segments3d(rbind(M1[dots$links[i, 1], ], M1[dots$links[i, 2], ]), col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
                    rgl::segments3d(rbind(M2[dots$links[i, 1], ], M2[dots$links[i, 2], ]), col = tarlinkcol[i], lty = tarlinklty[i], lwd = tarlinklwd[i])
                }
            }
        }
    }
}
