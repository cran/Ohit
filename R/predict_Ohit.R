#' Make predictions based on a fitted "Ohit" object
#'
#' @description
#' This function returns predictions from a fitted \code{"Ohit"} object.
#' @param object Fitted "Ohit" model object.
#' @param newX Matrix of new values for \code{X} at which predictions are to be made.
#' @return \item{pred_HDIC}{The predicted value based on the model determined by OGA+HDIC.}
#' @return \item{pred_Trim}{The predicted value based on the model determined by OGA+HDIC+Trim.}
#' @author Hai-Tang Chiou, Ching-Kang Ing and Tze Leung Lai.
#' @references Ing, C.-K. and Lai, T. L. (2011). A stepwise regression method and consistent model selection for high-dimensional sparse linear models. \emph{Statistica Sinica}, \strong{21}, 1473--1513.
#' @examples
#' # Example setup (Example 3 in Section 5 of Ing and Lai (2011))
#' n = 410
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
#' # with intercept
#' fit1 = Ohit(X[1:400, ], y[1:400])
#' predict_Ohit(fit1, rbind(X[401:401, ]))
#' predict_Ohit(fit1, X[401:410, ])
#' # without intercept
#' fit2 = Ohit(X[1:400, ], y[1:400], intercept = FALSE)
#' predict_Ohit(fit2, rbind(X[401:401, ]))
#' predict_Ohit(fit2, X[401:410, ])
#' @export

predict_Ohit <- function(object, newX){
  if (!is.matrix(newX)) stop("newX should be a matrix")
  if (ncol(newX) != object$p) stop(paste("the number of columns of newX is not equal to ", object$p, sep = ""))

  if (length(object$J_HDIC) == nrow(object$betahat_HDIC$coefficients)){
    pred_HDIC = newX[, object$J_HDIC] %*% as.matrix(object$betahat_HDIC$coefficients[, 1])
  }else if (nrow(newX) == 1){
    pred_HDIC = sum(c(1, newX[, object$J_HDIC]) * object$betahat_HDIC$coefficients[, 1])
  }else{
    pred_HDIC = cbind(rep(1, nrow(newX)), newX[, object$J_HDIC]) %*% as.matrix(object$betahat_HDIC$coefficients[, 1])
  }

  if (length(object$J_Trim) == nrow(object$betahat_Trim$coefficients)){
    pred_Trim = newX[, object$J_Trim] %*% as.matrix(object$betahat_Trim$coefficients[, 1])
  }else if (nrow(newX) == 1){
    pred_Trim = sum(c(1, newX[, object$J_Trim]) * object$betahat_Trim$coefficients[, 1])
  }else{
    pred_Trim = cbind(rep(1, nrow(newX)), newX[, object$J_Trim]) %*% as.matrix(object$betahat_Trim$coefficients[, 1])
  }

  return(list("pred_HDIC" = as.vector(pred_HDIC), "pred_Trim" = as.vector(pred_Trim)))
}

