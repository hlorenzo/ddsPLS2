#' Function to predict from ddsPLS objects
#'
#' @param x ddsPLS object
#' @param X_test matrix, a test data-set. If is "NULL", the default value, the predicted values for the train test are returned
#' @param ... Other parameters
#'
#' @export
#'
#' @seealso \code{\link{ddsPLS}}, \code{\link{plot.ddsPLS}}, \code{\link{summary.ddsPLS}}
#'
#' @useDynLib ddsPLS
predict.ddsPLS <- function(x,X_test=NULL,...){
  if(is.null(X_test)){
    y_est <- x$Y_est
  }else{
    n_test <- nrow(X_test)
    y_est <- (X_test-matrix(rep(x$model$muX,n_test),nrow = n_test,byrow = T))%*%x$model$B
    y_est <- y_est + matrix(rep(x$model$muY,n_test),nrow = n_test,byrow = T)
  }
  y_est
}
