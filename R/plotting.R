#' Screeplot of the Generalized Matrix Decomposition
#' @param fit An object of class "gmd" that is the output from \code{\link{GMD}}.
#' @details The function screeplot.gmd plots the proportion of the variances explained by each GMD component of an object of class "gmd".
#' @seealso \code{\link{GMD}}
#' @author Parker Knight and Yue Wang \email{ywang2310@fredhutch.org}
#' @references Allen, G. I., L. Grosenick, and J. Taylor (2014). A generalized least-square matrix decom- position. Journal of the American Statistical Association 109(505), 145–159.

#' @examples \dontrun{
#'    X = matrix(rnorm(1000), 50, 20)
#'    autocorr.mat <- function(p, rho) {
#'      mat <- diag(p)
#'      return(rho^abs(row(mat)-col(mat)))
#'      }
#'    H = autocorr.mat(50, 0.6)
#'    Q = autocorr.mat(20, 0.7)
#'    GMD.fit = GMD(X, H, Q, 2)
#'    screeplot(GMD.fit)
#'    }


#' @export
screeplot.gmd = function(fit){

  p = dim(fit$V.s)[1]

  D = fit$D
  X = fit$X
  H = fit$H
  Q = fit$Q
  K = length(D)

  D.plot = D^2/sum(diag(t(X)%*%H%*%X%*%Q))

  plot(1:K, D.plot[1:K], xaxt = 'n', ylab = 'Percentage of variance explained', type = 'b', xlab = ' ', ylim = c(0,ceiling(D.plot[1]*10)/10), cex.axis = 1.2, cex.lab = 1.2)
  axis(1,1:K, paste0("PC", 1:K))

}

#' The GMD-biplot
#' @description Biplots based on generelized matrix decomposition.
#' @param x An object of class "gmd" that is the output from \code{\link{GMD}}.
#' @param ... Optional arguments (see Details).
#' @details The optional arguments that can be passed in the biplot are
#' \itemize{
#'
#'  \item{\code{index}}{ a numeric vector specifying which variables(arrows) to display. The default is to plot the three longest arrows.}
#'  \item{\code{names}}{ a string vector specifying the names of the arrows. The default name of the i-th variable is Vi.}
#'  \item{\code{sample.col}}{ The color of the sample dots. The default is grey.}
#'  \item{\code{sample.pch}}{ Either an integer specifying a symbol or a single character to be used in plotting sample points. The default is 19.}
#'  \item{\code{arrow.col}}{ The color of the arrows. The default is orange.}
#'  \item{\code{arrow.cex}}{ A numeric value giving the size of the labels of the arrows. The default is 1.}
#'
#' }
#' @seealso \code{\link{GMD}}, \code{\link{screeplot.gmd}}
#' @author Parker Knight and Yue Wang \email{ywang2310@fredhutch.org}
#' @references Yue Wang, Timothy Randolph, Ali Shojaie and Jing Ma (2019). The GMD-biplot and its application to microbiome data.

#' @examples \dontrun{
#'    X = matrix(rnorm(1000), 50, 20)
#'    autocorr.mat <- function(p, rho) {
#'      mat <- diag(p)
#'      return(rho^abs(row(mat)-col(mat)))
#'      }
#'    H = autocorr.mat(50, 0.6)
#'    Q = autocorr.mat(20, 0.7)
#'    GMD.fit = GMD(X, H, Q, 2)
#'    screeplot(GMD.fit)
#'    biplot(GMD.fit)
#'    }


#' @export
biplot.gmd <- function(x, ...)
{
    biplot.gmd.body(x, ...)
}

biplot.gmd.body = function(fit, index, names, sample.col, sample.pch, arrow.col, arrow.cex){

  U = fit$U
  D = fit$D
  V = fit$V

  U = U[,order(D, decreasing = T)]
  V = V[,order(D, decreasing = T)]

  k1 = order(D, decreasing = T)[1]
  k2 = order(D, decreasing = T)[2]

  D = sort(D, decreasing = T)

  eta = U%*%diag(D)


  max.xlab = max(abs(eta[,1]))
  max.ylab = max(abs(eta[,2]))



  if(missing(sample.col)){sample.col = 'grey50'}
  if(missing(arrow.col)){arrow.col = 'orange'}
  if(missing(sample.pch)){sample.pch = 19}
  if(missing(arrow.cex)){arrow.cex = 1}

  plot(eta[,1], eta[,2], xlab = paste0('PC',k1), ylab = paste0('PC',k2), pch = sample.pch, xlim = c(-1.1*max.xlab, 1.1*max.xlab ), ylim =  c(-1.1*max.ylab, 1.1*max.ylab) , col = sample.col, cex.axis = 1.2, cex.lab = 1.2)
  xaxp = axTicks(1)
  yaxp = axTicks(2)


  # only plot these user-specified variables

  #calculate coordinates
  Q = fit$Q
  V.plot = Q%*%V
  arrow.x = V.plot[,1]
  arrow.y = V.plot[,2]

  # original code (equivalent!)
  #arrow.x = diag(rep(1, dim(V)[1]))%*%Q%*%V[,1]
  #arrow.y = diag(rep(1, dim(V)[1]))%*%Q%*%V[,2]

  if(missing(names)){names = paste0("V", index)}
  gmd.order = order(rowSums(V^2), decreasing = T)

  if(missing(index)){index = gmd.order[1:3]}

  max.xarrow = max(abs(arrow.x))
  max.yarrow = max(abs(arrow.y))
  xratio = max.xarrow/max.xlab
  yratio = max.yarrow/max.ylab


  xsci = as.numeric(unlist(strsplit(formatC(xratio, format = 'e'),"e")))
  xlab.arrow = round(xaxp*xsci[1]*10^(xsci[2]), digits = 2)


  ysci = as.numeric(unlist(strsplit(formatC(yratio, format = 'e'),"e")))
  ylab.arrow = round(yaxp*ysci[1]*10^(ysci[2]), digits = 2)

  iter = 1
  for(i in index){

    # if(arrow.x[i]^2 + arrow.y[i]^2 >= 0.1^2){

    arrows(x0 = 0,y0 = 0,x1 = arrow.x[i]/xratio, y1 = arrow.y[i]/yratio, length = 0.05, col = arrow.col)
    if(!missing(names)){
      text(arrow.x[i]/xratio, arrow.y[i]/yratio*1.1, names[iter], cex = 1, col = 'black')
    }

    #}

    iter = iter + 1
  }

  # add new axis
  axis(3, at = xaxp, labels = as.character(xlab.arrow), cex.axis = 1.2)
  axis(4, at = yaxp, labels = as.character(ylab.arrow), cex.axis = 1.2)

  points(eta[,1], eta[,2], col = sample.col, pch = sample.pch)
}
