#' Function to plot bootstrap performance results of the ddsPLS algorithm.
#'
#' @param x A ddsPLS object
#' @param type The type of graphics. One of "criterion" (default), "total",
#' "prop", "predict", "Q2r", "Q2", "R2r", "R2", "weightsX" or "weightsY".
#' @param digits double. Rounding of the written explained variance.
#' @param las integer. Rotation of label names.
#' @param cex.names double. Size factor for variable names.
#' @param ... Other plotting parameters to affect the plot.
#'
#' @importFrom graphics layout
#'
#' @export
#'
#' @useDynLib ddsPLS
plot.ddsPLS <- function(x,type="criterion",
                        digits=1,
                        legend.position="topright",
                        horiz=TRUE,
                        cex.names=1,mar=c(5, 4, 4, 2) + 0.1,
                        ...){
  ## Reset personnal plot par() settings
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  ## -----------------------------------
  h_opt <- x$R
  if(h_opt>0){
    q <- ncol(x$Y_obs)
    lambdas <- x$results$lambdas
    if(type=="total"){
      layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
    }
    if(h_opt>1){
      R2mean_diff_Q2mean <- matrix(do.call(cbind,x$results$R2mean_diff_Q2mean)[,1:h_opt],ncol = h_opt)
      Q2hmean <- matrix(do.call(cbind,x$results$Q2hmean)[,1:h_opt],ncol = h_opt)
      Q2mean <- matrix(do.call(cbind,x$results$Q2mean)[,1:h_opt],ncol = h_opt)
      R2hmean <- matrix(do.call(cbind,x$results$R2hmean)[,1:h_opt],ncol = h_opt)
      R2mean <- matrix(do.call(cbind,x$results$R2mean)[,1:h_opt],ncol = h_opt)
      PropQ2hPos <- matrix(do.call(cbind,x$results$PropQ2hPos)[,1:h_opt],ncol = h_opt)
    }else{
      R2mean_diff_Q2mean <- matrix(x$results$R2mean_diff_Q2mean[[1]],ncol = 1)
      Q2hmean <- matrix(x$results$Q2hmean[[1]],ncol = 1)
      Q2mean <- matrix(x$results$Q2mean[[1]],ncol = 1)
      R2hmean <- matrix(x$results$R2hmean[[1]],ncol = 1)
      R2mean <- matrix(x$results$R2mean[[1]],ncol = 1)
      PropQ2hPos <- matrix(x$results$PropQ2hPos[[1]],ncol = 1)
    }
    if(type %in% c("total","criterion")){
      # Plot of R2-Q2
      matplot(lambdas,R2mean_diff_Q2mean,type = "l",ylab="",xlab=expression(lambda),
              main=bquote(bar(R)[B]^2-bar(Q)[B]^2))
      for(s in 1:h_opt){
        points(lambdas,R2mean_diff_Q2mean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
        points(x$lambda[s],R2mean_diff_Q2mean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
      }
      legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
             col = 1:h_opt,pch=16,bty = "n",
             title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
    }
    if(type %in% c("total","prop")){
      # Plot of Prop of positive Q2h
      matplot(lambdas,PropQ2hPos,type = "l",ylab="",xlab=expression(lambda),
              main=bquote("Proportion of models with positive"~Q["b,r"]^2))
      abline(h=((1:10)/10)[-5],col="gray80",lwd=0.5,lty=3)
      abline(h=5/10,col="gray60",lwd=0.7,lty=1)
      text(min(lambdas),1/2,labels = "1/2",pos = 4,col="gray40")
      for(s in 1:h_opt){
        points(lambdas,PropQ2hPos[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
        points(x$lambda[s],PropQ2hPos[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
      }
      legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
             col = 1:h_opt,pch=16,bty = "n",
             title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
      q <- ncol(x$Y_obs)
    }
    if(type %in% c("total","predict")){
      # Predicted versus observed
      matplot(x$Y_obs,x$Y_est,pch=1:q,type="p",xlab="Observed",ylab="Predicted",
              main="Predicted versus observed")
      abline(0,1)
      if(is.null(colnames(x$Y_obs))){
        nono <- paste("Y",1:q," (",round(x$varExplained$PerY,digits = digits),"%)",sep="")
      } else {
        nono <- paste(colnames(x$Y_obs)," (",round(x$varExplained$PerY,digits = digits),"%)",sep="")
      }
      legend(legend.position,nono,col = 1:q,pch=1:q,bty = "n")
    }
    if(type %in% c("total","selection")){
      par(mfrow=c(2,1))
      selX <- x$Selection$X
      selY <- x$Selection$Y
      if(is.null(rownames(x$model$U))){
        names(selX) <- paste("X",1:p,sep="")
      }
      plot(selX,ylab="",
           main=paste("Which variables are selected in block X ?"),
           yaxt="n",xaxt="n",ylim=c(0,1),col=1+selX,pch=16+selX)
      axis(2,at = c(0,1),labels = c("No","Yes"),las=2)
      axis(1,at = 1:length(selX),labels = names(selX),
           las=las,cex.axis=cex.names)
      if(is.null(rownames(x$model$V))){
        names(selY) <- paste("Y",1:q,sep="")
      }
      plot(selY,ylab="",
           main=paste("Which variables are selected in block Y ?"),
           yaxt="n",xaxt="n",ylim=c(0,1),col=1+selY,pch=16+selY)
      axis(2,at = c(0,1),labels = c("No","Yes"),las=2)
      axis(1,at = 1:length(selY),labels = names(selY),
           las=las,cex.axis=cex.names)
    }
    if(type %in% c("total","Q2r")){
      # PLot of Q2_h
      matplot(lambdas,Q2hmean,type = "l",ylab="",xlab=expression(lambda),
              main=bquote(bar(Q)["B,r"]^2))
      for(s in 1:h_opt){  for(s in 1:h_opt){
        points(lambdas,Q2hmean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
        points(x$lambda[s],Q2hmean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
      }
        points(x$lambda[s],Q2hmean[which(lambdas==x$lambda[s]),s],pch=16,col=s)
      }
      abline(h=x$lowExplainedVariance,lwd=2,lty=2)
      legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
             col = 1:h_opt,pch=16,bty = "n",
             title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
    }
    if(type %in% c("R2r")){
      # PLot of R2_h
      matplot(lambdas,R2hmean,type = "l",ylab="",xlab=expression(lambda),
              main=bquote(bar(R)["B,r"]^2))
      for(s in 1:h_opt){  for(s in 1:h_opt){
        points(lambdas,R2hmean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
        points(x$lambda[s],R2hmean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
      }
        points(x$lambda[s],R2hmean[which(lambdas==x$lambda[s]),s],pch=16,col=s)
      }
      legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
             col = 1:h_opt,pch=16,bty = "n",
             title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
    }
    if(type %in% c("Q2")){
      # PLot of Q2
      matplot(lambdas,Q2mean,type = "l",ylab="",xlab=expression(lambda),
              main=bquote(bar(Q)["B"]^2))
      for(s in 1:h_opt){  for(s in 1:h_opt){
        points(lambdas,Q2mean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
        points(x$lambda[s],Q2mean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
      }
        points(x$lambda[s],Q2mean[which(lambdas==x$lambda[s]),s],pch=16,col=s)
      }
      abline(h=0,lwd=2,lty=2)
      legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
             col = 1:h_opt,pch=16,bty = "n",
             title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
    }
    if(type %in% c("R2")){
      # PLot of R2
      matplot(lambdas,R2mean,type = "l",ylab="",xlab=expression(lambda),
              main=bquote(bar(R)["B"]^2))
      for(s in 1:h_opt){  for(s in 1:h_opt){
        points(lambdas,R2mean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
        points(x$lambda[s],R2mean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
      }
        points(x$lambda[s],R2mean[which(lambdas==x$lambda[s]),s],pch=16,col=s)
      }
      abline(h=0,lwd=2,lty=2)
      legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
             col = 1:h_opt,pch=16,bty = "n",
             title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
    }
    if(type == "weightX" | type == "weightY"){
      q <- ncol(x$Y_obs)
      # layout(maLayout)
      par(mfrow=c(1,h_opt),mar=mar)
      if(is.null(colnames(x$Y_obs))){
        colnames(x$Y_obs) <- paste("Y",1:q,sep="")
      }
      if(is.null(rownames(x$model$U))){
        p <- nrow(x$model$U)
        rownames(x$model$U) <- paste("X",1:p,sep="")
      }
      for(s in 1:h_opt){
        if(type == "weightX"){
          popo <- t(x$model$U)[s,,drop=F]
          barplot(popo,
                  xlim=c(-1,1)*max(abs(popo)),
                  horiz = horiz,axis.lty = 1,las=2,
                  main=paste("X part, Comp.",s),
                  cex.names = cex.names)
        }else{
          popo <- t(x$model$V)[s,,drop=F]
          colnames(popo) <- paste(
            colnames(popo)," (",
            round(x$varExplained$PerYPerComp$Comp[s,],
                  digits = digits ),
            "%)",sep="")
          barplot(as.vector(popo),
                  xlim=c(-1,1)*max(abs(popo)),
                  cex.names = cex.names,horiz = horiz,
                  names=colnames(popo),col=(1:q),
                  axis.lty = 1,las=2,
                  main=paste("Y part, Comp.",s))
        }
        abline(v = c(-10:10)/10,lty=2,col="gray")
      }
    }
  }
}
