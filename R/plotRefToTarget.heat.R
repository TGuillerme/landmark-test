#' @title Plot shape differences as PlotRefToTarget with heat
#'
#' @description Modification of PlotRefToTarget to allow to add heat maps (only works for \code{methods = c("points", "vector")})
#'
#' @param M1 Matrix of landmark coordinates for the first (reference) specimen
#' @param M2 Matrix of landmark coordinates for the second (target) specimen
#' @param method Method used to visualize shape difference: either \code{"vector"} (default) or \code{"points"}. For \code{"surface"} or \code{"TPS"}, use \code{\link[geomorph]{plotRefToTarget}}.
#' @param plotRefToTarget.args additional arguments to be passed to \code{\link[geomorph]{plotRefTarget}}.
#' @param ... Additional parameters to \code{\link[graphics]{plot}} or \code{\link[rgl]{plot3d}}.
#' @param col Either a single color value (\code{"character"}), vector of values or function for colouring both points and vectors; or a list of two of any of these three elements. The first argument is passed to the points and the second to the vectors (\code{default = list("grey", "black")} for grey points and black vectors - see details).
#' @param col.val Optional, if \code{col} is a function for assigning colours (or a list of two functions), which values to pass to the function to assign the colors).
#' @param pt.size Optional, the size of the points (\code{default = 1}). This argument can overwrite \code{gridPars$pt.size} passed to \code{...}.
#' 
#' @details
#' \code{col} arguments are passed to \code{\link[graphics]{plot}} or \code{\link[rgl]{plot3d}} as follows:
#' \itemize{
#'    \item if col is a single value or a single vector, it passes it to both points and colors arguments;
#'    \item if col is a single function, it applies it to both points and colors arguments based on the argument \code{col.val};
#'    \item if col is a list of two values, vectors or functions, it will apply the first element to the points and the second to the vectors the same way as described above;
#' }
#' 
#' 
#' @examples
#' ## Loading the geomorph dataset
#' require(geomorph)
#' data(plethodon)
#' 
#' ## Performing the Procrustes superimposition
#' proc_super <- gpagen(plethodon$land, print.progress = FALSE)
#' 
#' ## Getting the range of variation
#' variation <- variation.range(proc_super, return.ID = TRUE)
#' 
#' ## Selecting the coordinates and the variation vector
#' M1 <- proc_super$coords[, , variation$min.max[1]]
#' M2 <- proc_super$coords[, , variation$min.max[2]]
#' var_val <- variation$range[, 1]
#'  
#' ## Plot the variation
#' plotRefToTarget.heat(M1, M2, method = "vector")
#' 
#' ## A weird looking plethodon
#' plotRefToTarget.heat(M1, M2, method = "points", col = list(rainbow, "pink"), col.val = var_val,
#'                 pt.size = 2.5, plotRefToTarget.args = list(mag = 3, outline = plethodon$outline))
#' 
#' \dontrun{
#' ## Loading the scallops 3D data from geomorph
#' require(geomorph)
#' data(scallops)
#' 
#' ## Procrustes superimposition
#' procrustes <- gpagen(scallops$coorddata)
#' 
#' ## Getting the range of variation
#' variation <- variation.range(procrustes, return.ID = TRUE)
#' 
#' ## Selecting the coordinates and the variation vector
#' M1 <- procrustes$coords[, , variation$min.max[1]]
#' M2 <- procrustes$coords[, , variation$min.max[2]]
#' var_val <- variation$range[, 1]
#'  
#' ## Plot the variation in 3D
#' plotRefToTarget.heat(M1, M2, method = "vector",  col.val = var_val,
#'                      col = list(grDevices::heat.colors, "grey"))
#' }
#' 
#' @seealso \code{\link{variation.range}}
#'  
#' @author Dean Adams, Emma Sherratt & Michael Collyer - modified by Thomas Guillerme
#' 
#' @export 
#' @importFrom graphics plot segments points text arrows
#' @importFrom rgl plot3d text3d

plotRefToTarget.heat <- function(M1, M2, method = "vector", plotRefToTarget.args = list(...), ..., col = list("grey", "black"), col.val, pt.size) {



  ## Get the methods (defaults)
  all_methods <- c("vector", "points")
  check.method(method, all_methods, "method")


  # method <- match.arg(method)

  ## Manage default arguments from plotRefToTarget
  dots <- plotRefToTarget.args

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
  if(is.null(dots$gridPars)) {
    gP = gridPar() 
  } else {
    gP = dots$gridPars
  }

  ## Check pt.size
  if(!missing(pt.size)) {
    gP$pt.size <- pt.size
  }

  k <- dim(M1)[2]

  dots$mag <- (dots$mag - 1)
  M2 <- M2 + (M2 - M1) * dots$mag
  
  limits = function(x, s) { 
      r = range(x)
      rc = scale(r, scale = FALSE)
      l = mean(r) + s * rc
  }


  ## Check the colour arguments
  ## Is it in the right format?
  if(class(col) == "list") {
    col_classes <- unlist(lapply(col, class))
      if(!any(col_classes %in% c("character", "function"))) {
        stop("col argument must be a single value, vector (of the same number of rows as M1) or function or a list of two of any of the former.")
      }
    col_list <- TRUE
    ## Should have two elements
    if(length(col) > 2) {
      stop("col argument must be a single value, vector (of the same number of rows as M1) or function or a list of two of any of the former.")
    }
  } else {
    if(!any(class(col) %in% c("character", "function"))) {
      stop("col argument must be a single value, vector (of the same number of rows as M1) or function or a list of two of any of the former.")
    }
    col_list <- FALSE
    ## Converting into a list for hanging
    col <- list(col)
  }

  ## Checking the arguments
  col_fun_check <- rep(FALSE, length(col))
  for(sub in 1:length(col)) {
    ## Get the class
    class_col <- class(col[[sub]])
    ## Check the argument
    if(class_col != "function") {
      if(class_col == "character") {
        if(length(col[[sub]]) > 1) {
          if(nrow(M1) < length(col[[sub]])) {
            stop("col argument must be a single value, vector (of the same number of rows as M1) or function or a list of two of any of the former.")
          }
        } else {
          col[[sub]] <- rep(col[[sub]], nrow(M1))
        }
      } else {
        stop("col argument must be a single value, vector (of the same number of rows as M1) or function or a list of two of any of the former.")
      }
    } else {
      ## Colour is a function
      col_fun_check[sub] <- TRUE
    }
  }

  ## Checking the col.val argument
  if(!missing(col.val)) {
    if(class(col.val) == "list") {
      if(length(col.val) > 2) {
        stop("col.val should be a list of two vectors or a vector, all of the same length as the number of rows in M1 and M2.")
      }
    } else {
      col.val <- list(col.val, col.val)
    }

    ## Check the length
    col_val_length <- unlist(lapply(col.val, length))
    if(any(col_val_length != nrow(M1))) {
      stop("col.val should be a list of two vectors or a vector, all of the same length as the number of rows in M1 and M2.")
    } else {
      val.points <- col.val[[1]]
      val.vector <- col.val[[2]]
    }
  } else {
    if(any(col_fun_check)) {
      stop("col.val argument is missing with no default.")
    }
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

  ## Handling the colour arguments
  if(col_list) {
    if(col_fun_check[1]) {
      ## Apply the function
      col_points <- get.colour.gradient(val.points, col[[1]])
    } else {
      col_points <- col[[1]]
    }
    if(col_fun_check[2]) {
      ## Apply the function
      col_vector <- get.colour.gradient(val.vector, col[[2]])
    } else {
      col_vector <- col[[2]]
    }
  } else {
    if(any(col_fun_check)) {
      col_points <- col_vector <- get.colour.gradient(val.points, col[[1]])
    } else {
      col_points <- col_vector <- col[[1]]
    }
  }

  ## 2D plots
  if(k == 2) {

    if(method == "vector") {
        if(dots$axes == TRUE) {
          graphics::plot(M1, asp = 1, type = "n", xlab = "x", ylab = "y", xlim = limits(M1[, 1], 1.25), ylim = limits(M1[, 2], 1.25), ...)
        } else {
          graphics::plot(M1, asp = 1, type = "n", xlab = "", ylab = "", xlim = limits(M1[, 1], 1.25), axes = FALSE, ylim = limits(M1[, 2], 1.25), ...)
        }
        if(is.null(dots$links) == FALSE) {
            linkcol <- rep(gP$link.col, nrow(dots$links))[1:nrow(dots$links)]
            linklwd <- rep(gP$link.lwd, nrow(dots$links))[1:nrow(dots$links)]
            linklty <- rep(gP$link.lty, nrow(dots$links))[1:nrow(dots$links)]
            
            for (i in 1:nrow(dots$links)) {
                graphics::segments(M2[dots$links[i, 1], 1], M2[dots$links[i, 1], 2], M2[dots$links[i, 2], 1], M2[dots$links[i, 2], 2], col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
            }
        }
        if(dots$label == TRUE) {
          graphics::text(M1, label = paste(1:dim(M1)[1]), adj = gP$txt.adj, pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)
        }

        graphics::arrows(M1[, 1], M1[, 2], M2[, 1], M2[, 2], length = 0.075, lwd = 2, col = col_vector)

        graphics::points(M1, pch = 21, bg = col_points, cex = gP$pt.size, col = gP$pt.bg)
    }

    if(method == "points") {
        if(dots$axes == TRUE) {
          graphics::plot(M1, asp = 1, pch = 21, type = "n", xlim = limits(M1[, 1], 1.25), ylim = limits(M1[, 2], 1.25), xlab = "x", ylab = "y", ...)
        }
        if(dots$axes == FALSE) {
            graphics::plot(M1, asp = 1, pch = 21, type = "n", xlim = limits(M1[, 1], 1.25), axes = FALSE, ylim = limits(M1[, 2], 1.25), xlab = "", ylab = "", ...)
        }
        if(dots$label == TRUE) {text(M1, label = paste(1:dim(M1)[1]), adj = gP$txt.adj, pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)}
        if(!is.null(dots$outline)) {

          ## Internal function from geomorph
          tps2d<-function(M, matr, matt)
          {p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
          P<-matrix(NA, p, p)
          for (i in 1:p)
          {for (j in 1:p){
            r2<-sum((matr[i,]-matr[j,])^2)
            P[i,j]<- r2*log(r2)}}
          P[which(is.na(P))]<-0
          Q<-cbind(1, matr)
          L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
          m2<-rbind(matt, matrix(0, 3, 2))
          fast.ginv <- function(X, tol = sqrt(.Machine$double.eps)){
            k <- ncol(X)
            Xsvd <- La.svd(X, k, k)
            Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
            rtu <-((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
            v <-t(Xsvd$vt)[, Positive, drop = FALSE]
            v%*%rtu
          }
          fast.solve <- function(x) if(det(x) > 1e-8) qr.solve(x) else fast.ginv(x)
          coefx<-fast.solve(L)%*%m2[,1]
          coefy<-fast.solve(L)%*%m2[,2]
          fx<-function(matr, M, coef)
          {Xn<-numeric(q)
          for (i in 1:q)
          {Z<-apply((matr-matrix(M[i,],p,2,byrow=TRUE))^2,1,sum)
          Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))}
          Xn}
          matg<-matrix(NA, q, 2)
          matg[,1]<-fx(matr, M, coefx)
          matg[,2]<-fx(matr, M, coefy)
          matg}

            curve.warp <- tps2d(dots$outline, M1, M2)
            graphics::points(dots$outline, pch = 19, cex = gP$out.cex, col = gP$out.col) 
            graphics::points(curve.warp, pch = 19, cex = gP$tar.out.cex, col = gP$tar.out.col) 
        }
        if(is.null(dots$links) == FALSE) {
            linkcol <- rep(gP$link.col, nrow(dots$links))[1:nrow(dots$links)]
            linklwd <- rep(gP$link.lwd, nrow(dots$links))[1:nrow(dots$links)]
            linklty <- rep(gP$link.lty, nrow(dots$links))[1:nrow(dots$links)]
            tarlinkcol <- rep(gP$tar.link.col, nrow(dots$links))[1:nrow(dots$links)]
            tarlinklwd <- rep(gP$tar.link.lwd, nrow(dots$links))[1:nrow(dots$links)]
            tarlinklty <- rep(gP$tar.link.lty, nrow(dots$links))[1:nrow(dots$links)]
            for (i in 1:nrow(dots$links)) {
                graphics::segments(M1[dots$links[i, 1], 1], M1[dots$links[i, 1], 2], M1[dots$links[i, 2], 1], M1[dots$links[i, 2], 2], col = linkcol[i],lty = linklty[i], lwd = linklwd[i])
                graphics::segments(M2[dots$links[i, 1], 1], M2[dots$links[i, 1], 2], M2[dots$links[i, 2], 1], M2[dots$links[i, 2], 2], col = tarlinkcol[i], lty = tarlinklty[i], lwd = tarlinklwd[i])
            }
        }
        graphics::points(M2, pch = 21, bg = col_points, cex = gP$tar.pt.size)
        graphics::points(M1, pch = 21, bg = col_points, cex = gP$pt.size)
      }

      return(invisible())

    }

    ## 3D plots
    if(k == 3) {

        if(method == "vector") {

            if(dots$axes == TRUE) {
              rgl::plot3d(M1, type = "s", col = col_points, size = gP$pt.size, aspect = FALSE, ...)
            }
            if(dots$axes == FALSE) {
                rgl::plot3d(M1, type = "s", col = col_points, size = gP$pt.size, aspect = FALSE, xlab = "", ylab = "", zlab = "", axes = FALSE, ...)
            }
            if(dots$label == TRUE) {
              rgl::text3d(M1, texts = paste(1:dim(M1)[1]), adj = (gP$txt.adj + gP$pt.size), pos = (gP$txt.pos + gP$pt.size), cex = gP$txt.cex, col = gP$txt.col)
            } 

            ## Function for adding the segments
            add.segment <- function(i, M1, M2, lwd = 2, col) {
              segments3d(rbind(M1[i, ], M2[i, ]), lwd = lwd, col = col[i])
            }

            silent <- sapply(1:nrow(M1), add.segment, M1 = M1, M2 = M2, col = col_vector)

            # for (i in 1:nrow(M1)) {
            #   silent <- rgl::segments3d(rbind(M1[i, ], M2[i, ]), lwd = 2, col = col_vector[i])
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
              rgl::plot3d(M2, type = "s", col = col_vector, size = gP$tar.pt.size, add = TRUE)
            }
            if(dots$axes == FALSE) {
                rgl::plot3d(M1, type = "s", col = col_points, size = gP$pt.size, aspect = FALSE, xlab = "", ylab = "", zlab = "", axes = F, ...)
                rgl::plot3d(M2, type = "s", col = col_vector, size = gP$tar.pt.size, add = TRUE)
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
                    segments3d(rbind(M1[dots$links[i, 1], ], M1[dots$links[i, 2], ]), col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
                    rgl::segments3d(rbind(M2[dots$links[i, 1], ], M2[dots$links[i, 2], ]), col = tarlinkcol[i], lty = tarlinklty[i], lwd = tarlinklwd[i])
                }
            }
        }
    
      return(invisible())

    }
}
