get_uv = function(X, H, Q, u_0, v_0){

  u = X%*%Q%*%v_0/as.numeric(sqrt(t(v_0)%*%t(Q)%*%t(X)%*%H%*%X%*%Q%*%v_0))
  v = t(X)%*%H%*%u/as.numeric(sqrt(t(u)%*%t(H)%*%X%*%Q%*%t(X)%*%H%*%u))

  return(list(u = u, v = v))

}

#' Generalized Matrix Decomposition
#'
#' Computes the generalized matrix decomposition of X.
#'
#' @param X An n x p data matrix.
#' @param H An n x n positive semi-definite similarity kernel.
#' @param Q An p x p positive semi-definite similarity kernel.
#' @param K a scalar specifying the dimension of GMD components (see Details).
#'
#' @details The K-dimensional GMD of X with respect to H and Q is given by X = USV^T, where
#' \deqn{(U, S, V) = argmin_{U,S,V} ||X - USV^T||^2_{H, Q},}
#' subject to \eqn{U^THU = I_K, V^TQV = I_K} and \eqn{diag(S) \ge 0}. Here, for any n x p matrix A,
#' ||A||^2_{H,Q} = tr(A^THAQ).
#'
#'
#'
#
#' @return A list of the GMD components U, S, V, H and Q (see Details).
#' @author Parker Knight and Yue Wang \email{ywang2310@fredhutch.org}
#' @references Allen, G. I., L. Grosenick, and J. Taylor (2014). A generalized least-square matrix decom- position. Journal of the American Statistical Association 109(505), 145–159.
#' @examples
#' \dontrun{
#'    X = matrix(rnorm(1000), 50, 20)
#'    autocorr.mat <- function(p, rho) {
#'      mat <- diag(p)
#'      return(rho^abs(row(mat)-col(mat)))
#'      }
#'    H = autocorr.mat(50, 0.6)
#'    Q = autocorr.mat(20, 0.7)
#'    GMD.fit = GMD(X, H, Q, 2)
#'    }
#'
#' @export
GMD =function(X, H, Q, K){

  n = dim(X)[1]
  p = dim(X)[2]
  # output matrix/vec
  U = matrix(0, n, K)
  V = matrix(0, p, K)
  D = rep(0,K)

  X_0 = X
  u_0 = c(1, rep(0,n-1))
  v_0 = c(1, rep(0,p-1))

  for(iter in 1:K){


    error = 1

    while(error > 1e-5){

      temp.uv = get_uv(X_0, H, Q, u_0, v_0)

      u = temp.uv$u
      v = temp.uv$v

      error = norm(u - u_0, "2") + norm(v - v_0, "2")
      #print(error)

      u_0 = u
      v_0 = v

    }

    U[,iter] = u
    V[,iter] = v
    d = t(u)%*%H%*%X_0%*%Q%*%v
    D[iter] =  d

    X_0 = X_0 - u%*%d%*%t(v)

  }

  gmd.output <- list(U = U, V = V, D = D, H = H, Q = Q, X = X)
  class(gmd.output) <- "gmd"
  return(gmd.output)

}
