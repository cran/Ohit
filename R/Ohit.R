#' Fit a high-dimensional linear regression model via OGA+HDIC+Trim
#'
#' @description
#' The first step is to sequentially select input variables via orthogonal greedy algorithm (OGA). The second step is to determine the number of OGA iterations using high-dimensional information criterion (HDIC). The third step is to remove irrelevant variables remaining in the second step using HDIC.
#' @import stats
#' @param X Input matrix of \code{n} rows and \code{p} columns.
#' @param y Response vector of length \code{n}.
#' @param Kn The number of OGA iterations. \code{Kn} must be a positive integer between \code{1} and \code{p}. Default is \code{Kn=max(1, min(floor(c1}*\code{sqrt(n/log(p))), p))}, where \code{c1} is a tuning parameter.
#' @param c1 The tuning parameter for the number of OGA iterations. Default is \code{c1=5}.
#' @param HDIC_Type High-dimensional information criterion. The value must be \code{"HDAIC"}, \code{"HDBIC"} or \code{"HDHQ"}. The formula is \code{n}*\code{log(rmse)+k_use}*\code{omega_n}*\code{log(p)} where \code{rmse} is the residual mean squared error and \code{k_use} is the number of variables used to fit the model. For \code{HDIC_Type="HDAIC"}, it is HDIC with \code{omega_n=c2}. For \code{HDIC_Type="HDBIC"}, it is HDIC with \code{omega_n=log(n)}. For \code{HDIC_Type="HDHQ"}, it is HDIC with \code{omega_n=c3}*\code{log(log(n))}. Default is \code{HDIC_Type="HDBIC"}.
#' @param c2 The tuning parameter for \code{HDIC_Type="HDAIC"}. Default is \code{c2=2}.
#' @param c3 The tuning parameter for \code{HDIC_Type="HDHQ"}. Default is \code{c3=2.01}.
#' @param intercept Should an intercept be fitted? Default is \code{intercept=TRUE}.
#' @return \item{n}{The number of observations.}
#' @return \item{p}{The number of input variables.}
#' @return \item{Kn}{The number of OGA iterations.}
#' @return \item{J_OGA}{The index set of Kn variables sequencially selected by OGA.}
#' @return \item{HDIC}{The HDIC values along the OGA path.}
#' @return \item{J_HDIC}{The index set of valuables determined by OGA+HDIC.}
#' @return \item{J_Trim}{The index set of valuables determined by OGA+HDIC+Trim.}
#' @return \item{betahat_HDIC}{The estimated regression coefficients of the model determined by OGA+HDIC.}
#' @return \item{betahat_Trim}{The estimated regression coefficients of the model determined by OGA+HDIC+Trim.}
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
#' # Fit a high-dimensional linear regression model via OGA+HDIC+Trim
#' Ohit(X, y, intercept = FALSE)
#' @export

Ohit <- function(X, y, Kn = NULL, c1 = 5, HDIC_Type = "HDBIC", c2 = 2, c3 = 2.01, intercept = TRUE){
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

  if ((HDIC_Type != "HDAIC") & (HDIC_Type != "HDBIC") & (HDIC_Type != "HDHQ")) stop("HDIC_Type should be \"HDAIC\", \"HDBIC\" or \"HDHQ\"")
  if (HDIC_Type == "HDAIC") omega_n = c2
  if (HDIC_Type == "HDBIC") omega_n = log(n)
  if (HDIC_Type == "HDHQ") omega_n = c3 * log(log(n))

  hdic = (n * log(sigma2hat))+((1:K) * omega_n * (log(p)))
  kn_hat = which.min(hdic)
  benchmark = hdic[kn_hat]
  J_HDIC = sort(Jhat[1:kn_hat])

  J_Trim = Jhat[1:kn_hat]
  trim_pos = rep(0, kn_hat)
  if (kn_hat > 1){
    for (l in 1:(kn_hat-1)){
      JDrop1 = J_Trim[-l]
      fit = lm(dy~.-1, data = data.frame(dX[, JDrop1]))
      uDrop1 = fit$residuals
      HDICDrop1 = (n * log(mean(uDrop1 ^ 2)))+((kn_hat - 1) * omega_n * (log(p)))
      if (HDICDrop1 > benchmark) trim_pos[l] = 1
    }
    trim_pos[kn_hat] = 1
    J_Trim = J_Trim[which(trim_pos==1)]
  }
  J_Trim = sort(J_Trim)

  X_HDIC = as.data.frame(as.matrix(X[, J_HDIC]))
  X_Trim = as.data.frame(as.matrix(X[, J_Trim]))
  X = data.frame(X)
  colnames(X_HDIC) = names(X)[J_HDIC]
  colnames(X_Trim) = names(X)[J_Trim]

  if (intercept == TRUE){
    fit_HDIC = lm(y~., data = X_HDIC)
    fit_Trim = lm(y~., data = X_Trim)
  }else{
    fit_HDIC = lm(y~.-1, data = X_HDIC)
    fit_Trim = lm(y~.-1, data = X_Trim)
  }

  betahat_HDIC = summary(fit_HDIC)
  betahat_Trim = summary(fit_Trim)

  return(list("n" = n, "p" = p, "Kn" = K, "J_OGA" = Jhat, "HDIC" = hdic, "J_HDIC" = J_HDIC, "J_Trim" = J_Trim, "betahat_HDIC" = betahat_HDIC, "betahat_Trim" = betahat_Trim))
}
