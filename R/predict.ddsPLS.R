#' Function to predict from ddsPLS objects
#'
#' @param x A ddsPLS object
#' @param X_test Optimal number of components to be plotted
#' @param ... Other parameters.
#'
#' @export
#'
#' @useDynLib ddsPLS
predict.ddsPLS <- function(x,X_test,...){
  n_test <- nrow(X_test)
  y_est <- (X_test-matrix(rep(x$model$muX,n_test),nrow = n_test,byrow = T))%*%x$model$B
  y_est + matrix(rep(x$model$muY,n_test),nrow = n_test,byrow = T)
}
