#' Orthogonal greedy algorithm
#'
#' @description
#' Select valuables via orthogonal greedy algorithm (OGA).
#' @import stats
#' @param X Input matrix of \code{n} rows and \code{p} columns.
#' @param y Response vector of length \code{n}.
#' @param Kn The number of OGA iterations. \code{Kn} must be a positive integer between \code{1} and \code{p}. Default is \code{Kn=max(1, min(floor(c1}*\code{sqrt(n/log(p))), p))}, where \code{c1} is a tuning parameter.
#' @param c1 The tuning parameter for the number of OGA iterations. Default is \code{c1=5}.
#' @return \item{n}{The number of observations.}
#' @return \item{p}{The number of input variables.}
#' @return \item{Kn}{The number of OGA iterations.}
#' @return \item{J_OGA}{The index set of \code{Kn} variables sequencially selected by OGA.}
#' @author Hai-Tang Chiou, Ching-Kang Ing and Tze Leung Lai.
#' @references Ing, C.-K. and Lai, T. L. (2011). A stepwise regression method and consistent model selection for high-dimensional sparse linear models. \emph{Statistica Sinica}, \strong{21}, 1473--1513.
#' @examples
#' # Example setup (Example 3 in Section 5 of Ing and Lai (2011))
#' n = 400
#' p = 4000
#' q = 10
#' beta_1q = c(3, 3.75, 4.5, 5.25, 6, 6.75, 7.5, 8.25, 9, 9.75)
#' b = sqrt(3/(4 * q))
#'
#' x_relevant = matrix(rnorm(n * q), n, q)
#' d = matrix(rnorm(n * (p - q), 0, 0.5), n, p - q)
#' x_relevant_sum = apply(x_relevant, 1, sum)
#' x_irrelevant = apply(d, 2, function(a) a + b * x_relevant_sum)
#' X = cbind(x_relevant, x_irrelevant)
#' epsilon = rnorm(n)
#' y = as.vector((x_relevant %*% beta_1q) + epsilon)
#'
#' # Select valuables via OGA
#' OGA(X, y)
#' @export

OGA <- function(X, y, Kn = NULL, c1 = 5){
  if (!is.vector(y)) stop("y should be a vector")
  if (!is.matrix(X)) stop("X should be a matrix")

  n = nrow(X)
  p = ncol(X)
  if (n != length(y)) stop("the number of observations in y is not equal to the number of rows of X")
  if (n == 1) stop("the sample size should be greater than 1")
  if (is.null(Kn)) K = max(1, min(floor(c1 * sqrt(n / log(p))), p))
  else{
    if ((Kn < 1) | (Kn > p)) stop(paste("Kn should between 1 and ", p, sep = ""))
    if ((Kn - floor(Kn)) != 0) stop("Kn should be a positive integer")
    K = Kn
  }

  dy = y - mean(y)
  dX = apply(X, 2, function(x) x - mean(x))

  Jhat = sigma2hat = rep(0, K)
  XJhat = matrix(0, n, K)
  u = as.matrix(dy)
  xnorms = sqrt(colSums((dX) ^ 2))

  aSSE = (abs(t(u) %*% dX) / xnorms)
  Jhat[1] = which.max(aSSE)
  XJhat[, 1] = (dX[, Jhat[1]] / sqrt(sum((dX[, Jhat[1]]) ^ 2)))
  u = u - XJhat[, 1] %*% t(XJhat[, 1]) %*% u
  sigma2hat[1] = mean(u ^ 2)

  if (K > 1){
    for (k in 2:K) {
      aSSE = (abs(t(u) %*% dX) / xnorms)
      aSSE[Jhat[1:(k-1)]] = 0
      Jhat[k] = which.max(aSSE)
      rq = dX[, Jhat[k]] - XJhat[, 1:(k-1)] %*% t(XJhat[, 1:(k-1)]) %*% dX[, Jhat[k]]
      XJhat[, k] = (rq / sqrt(sum((rq) ^ 2)))
      u = u - XJhat[, k] %*% t(XJhat[, k]) %*% u
      sigma2hat[k] = mean(u ^ 2)
    }
  }

  return(list("n" = n, "p" = p, "Kn" = K, "J_OGA" = Jhat))
}
