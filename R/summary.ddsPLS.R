#' Function to sum up bootstrap performance results of the ddsPLS algorithm.
#'
#' @param x A ddsPLS object.
#' @param return Wether or not to return the printed values, default to FALSE.
#' @param digits integer indicating the number of decimal places (round) to be used.
#' @param ... Other parameters to be taken into account.
#'
#' @export
#'
#' @useDynLib ddsPLS
summary.ddsPLS <- function(x,return=FALSE,
                           plotSelection=FALSE,las=1,cex.names=1,
                           digits=2,...){
  cat("                      ______________\n");
  cat("                     |    ddsPLS    |\n");
  cat("=====================----------------=====================\n");
  h_opt <- x$R
  if(h_opt>0){
    q <- ncol(x$Y_obs)
    p <- nrow(x$model$U)
    ## R2 and Q2 final values
    R2Q2 <- matrix(0,h_opt,5)
    for(h in 1:h_opt){
      best <- x$lambda[h]
      bestID <- which(x$results$lambdas==best)
      R2Q2[h,] <- c(
        best,
        x$results$R2mean[[h]][bestID],
        x$results$R2hmean[[h]][bestID],
        x$results$Q2mean[[h]][bestID],
        x$results$Q2hmean[[h]][bestID]
      )
    }
    colnames(R2Q2) <- c("lambda","R2","R2_r","Q2","Q2_r")
    rownames(R2Q2) <- paste("Comp.",1:h_opt)
    ## Explained variances for the components
    VAR <- do.call(rbind,x$varExplained[1:2])
    rownames(VAR) <- c("Per component ", "Cumulative ")
    colnames(VAR) <- paste("Comp.",1:h_opt)
    ## Explained variances for the response variables
    VARY <- matrix(x$varExplained$PerY,nrow=1)
    if(is.null(colnames(x$Y_obs))){
      colnames(VARY) <- paste("Y",1:q,sep="")
    } else {
      colnames(VARY) <- colnames(x$Y_obs)
    }
    rownames(VARY) <- ""
    ## Explained variances per response variables per component marginally
    VARYCompMarg <- matrix(x$varExplained$PerYPerComp$Comp,ncol=q)
    if(is.null(colnames(x$Y_obs))){
      colnames(VARYCompMarg) <- paste("Y",1:q,sep="")
    } else {
      colnames(VARYCompMarg) <- colnames(x$Y_obs)
    }
    rownames(VARYCompMarg) <- paste("Comp.",1:h_opt)
    ## ... in total
    VARYCompTotal <- matrix(x$varExplained$PerYPerComp$Cum,ncol=q)
    if(is.null(colnames(x$Y_obs))){
      colnames(VARYCompTotal) <- paste("Y",1:q,sep="")
    } else {
      colnames(VARYCompTotal) <- colnames(x$Y_obs)
    }
    rownames(VARYCompTotal) <- paste("Comp.",1:h_opt)
    ######
    ######
    cat(paste("The optimal ddsPLS model is built on",h_opt,"components\n\n"))
    ######
    cat(paste("The bootstrap quality statistics:\n",
              "---------------------------------\n",sep=""))
    print(round(R2Q2,digits))
    ######
    cat(paste("\n\nThe explained variance (in %):\n",
                "-----------------------\n",sep=""))
    cat(paste("\nIn total: ",round(x$varExplained$Cumu[h_opt],digits),
              "\n-  -  -  \n",sep=""))

    cat(paste("\nPer component or cumulated:\n",
              "-  -  -  -  -  -  -  -  -  \n",sep=""))
    print(round(VAR,digits))
    ######
    cat(paste("\nPer response variable:\n",
                "-  -  -  -  -  -  -  -\n",sep=""))
    print(round(VARY,digits))
    ######
    cat(paste("\nPer response variable per component:\n",
                "-  -  -  -  -  -  -  -  -  -  -  -  \n",
              sep=""))
    print(round(VARYCompMarg,digits))
    ######
    cat(paste("\n...and cumulated to:\n",
                "-  -  -  -  -  -  - \n",
              sep=""))
    print(round(VARYCompTotal,digits))
  }else{
    cat(paste("The optimal ddsPLS model is empty.\n"))
  }
  cat("=====================                =====================\n");
  cat("                     ================\n");
  if(return & h_opt>0){
    out <- list(R2Q2=R2Q2,VAR=VAR,VARY=VARY,
                VARYperComp=list(marginal=VARYCompMarg,total=VARYCompTotal))
    return(out)
  }
}