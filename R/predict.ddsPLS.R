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
predict.ddsPLS <- function(x,X_test=NULL,doDiagnosis=T,...){
  getDiagnoses <- function(x,X_test,y_test_est){
    n_test <- nrow(X_test)
    n_train <- nrow(x$X)
    X_test_center <- X_test-matrix(rep(x$model$muX,n_test),nrow = n_test,byrow = T)
    X_train_center <- x$X-matrix(rep(x$model$muX,n_train),nrow = n_train,byrow = T)
    ## On the predictions
    y_train_est <- x$Y_est
    sigma_y_inv <- unlist(lapply(x$model$sdY,function(ss){
      if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    y_mean_train <- matrix(rep(colMeans(y_train_est),n_train),nrow = n_train,byrow = T)
    y_mean_test <- matrix(rep(colMeans(y_train_est),n_test),nrow = n_test,byrow = T)
    dist_y_Mahalanobis_train <- (rowMeans(
      (y_train_est-y_mean_train)^2*matrix(rep(sigma_y_inv,n_train),nrow = n_train,byrow = T)^2))
    dist_y_Mahalanobis_test <- (rowMeans(
      (y_test_est-y_mean_test)^2*matrix(rep(sigma_y_inv,n_test),nrow = n_test,byrow = T)^2))
    ## On the scores
    t_test <- X_test_center%*%x$model$U_star
    t_train <- X_train_center%*%x$model$U_star
    sigma_t_inv <- unlist(lapply(sqrt(colMeans(t_train^2)),function(ss){
      if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    dist_score_Mahalanobis_train <- (rowMeans(
      t_train^2*matrix(rep(sigma_t_inv,n_train),nrow = n_train,byrow = T)^2))
    dist_score_Mahalanobis_test <- (rowMeans(
      t_test^2*matrix(rep(sigma_t_inv,n_test),nrow = n_test,byrow = T)^2))
    ## On the reconstruction of X
    X_test_center_est <- tcrossprod(t_test,x$model$P)
    X_train_center_est <- tcrossprod(t_train,x$model$P)
    epsilon_X_test <- X_test_center_est-X_test_center
    epsilon_X_train <- X_train_center_est-X_train_center
    # sigmaInv <- unlist(lapply(x$model$sdX,function(ss){
    #   if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    sigma_x_Inv <- unlist(lapply(sqrt(colMeans(epsilon_X_train^2)),function(ss){
      if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    dist_Mahalanobis_train <- (rowMeans(
      epsilon_X_train^2*matrix(rep(sigma_x_Inv,n_train),nrow = n_train,byrow = T)^2))
    dist_Mahalanobis_test <- (rowMeans(
      epsilon_X_test^2*matrix(rep(sigma_x_Inv,n_test),nrow = n_test,byrow = T)^2))
    list(# y=list(train=dist_y_Mahalanobis_train,test=dist_y_Mahalanobis_test),
         t=list(train=dist_score_Mahalanobis_train,test=dist_score_Mahalanobis_test),
         x=list(train=dist_Mahalanobis_train,test=dist_Mahalanobis_test))
  }
  diagnoses <- NULL
  if(is.null(X_test)){
    y_est <- x$Y_est
  }else{
    n_test <- nrow(X_test)
    y_est <- (X_test-matrix(rep(x$model$muX,n_test),nrow = n_test,byrow = T))%*%x$model$B
    y_est <- y_est + matrix(rep(x$model$muY,n_test),nrow = n_test,byrow = T)
    if(doDiagnosis){
      diagnoses <- getDiagnoses(x,X_test,y_est)
      plot(diagnoses$t$train,diagnoses$x$train,col="gray",pch=16,
           cex=1/2,xlab="t",ylab=expression(hat(x)))
      points(diagnoses$t$test,diagnoses$x$test,col="red",pch=16,cex=1)
    }
  }
  list(y_est=y_est,
       diagnoses=diagnoses)
}
