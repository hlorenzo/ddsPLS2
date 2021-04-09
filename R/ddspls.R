#' bootstrapWrap
#'
#' @param U
#' @param V
#' @param X
#' @param Y
#' @param lambdas
#' @param lambda_prev
#' @param R
#' @param n_B
#' @param doBoot
#' @param n
#' @param p
#' @param q
#' @param N_lambdas
#'
#' @return List
bootstrapWrap <- function(U,V,X,Y,lambdas,lambda_prev,
                          R,n_B,doBoot=TRUE,n,p,q,N_lambdas){
  res <- bootstrap_Rcpp(U,V,X,Y,lambdas,lambda_prev,
                        R,n_B,doBoot,n,p,q,N_lambdas)
  res
}

#' ddsPLS
#'
#' @param X
#' @param Y
#' @param lambdas
#' @param n_B
#' @param minBootProp
#' @param lowExplainedVariance
#' @param NCORES
#' @param errorMin
#' @param verbose
#'
#' @return
#' @export list
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @useDynLib ddsPLS
ddsPLS <- function(X,Y,lambdas,n_B,
                   minBootProp=0.0,
                   lowExplainedVariance=0.0,NCORES=1,errorMin=1e-9,verbose=FALSE){
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y)
  N_lambdas <- length(lambdas)
  # Standardize X and Y train and test.
  sdY <- apply(Y,2,sd)
  muY <- apply(Y,2,mean)
  Y_init = scale(Y);
  Y_init[which(is.na(Y_init))] <- 0
  RSS0 <- sum(scale(Y,scale = F)^2)
  RSS0_y <- apply(Y,2,function(yy){sum(scale(yy,scale = F)^2)})
  sdX <- apply(X,2,sd)
  muX <- apply(X,2,mean)
  X_init = scale(X);
  X_init[which(is.na(X_init))] <- 0
  sd_y_x_inv <- matrix(0,p,q)
  for(j in 1:q){
    sd_y_x_inv[,j] <- sdY[j]
  }
  for(i in 1:p){
    if(sdX[i]!=0){
      sd_y_x_inv[i,] <- 1/sdX[i]
    }else{
      sd_y_x_inv[i,] <- 0
    }
  }
  COVInit = crossprod(Y_init,X_init)/(n-1);
  maxCOVInit = max(abs(COVInit))
  lambda_prev <- rep(0,n)
  TEST <- rep(0,N_lambdas)
  test_lambdas <- list()
  nb_ValsOk <- 0
  U_out <- matrix(0,p,n); V0 <- matrix(0,q,n)
  B_previous <- matrix(0,p,q)
  R <- 1
  test <- TRUE
  R2Best <- rep(0,n)
  R2hBest <- rep(0,n)
  Q2Best <- rep(0,n)
  Q2hBest <- rep(0,n)
  explainedVar <- rep(0,n)
  Results <- list()
  varExplained = varExplainedTot = varExplained_y = varExplainedTot_y <- NULL
  Results$R2 = Results$R2h = Results$Q2 = Results$Q2h = Results$PropQ2hPos  <- list()
  Results$R2mean = Results$R2hmean = Results$Q2mean = Results$Q2hmean = Results$R2mean_diff_Q2mean <- list()
  h <- 0 ; bestID <- 0;
  Q2_previous <- -1e9 ; bestVal <- -1e9
  if (verbose) {
    cat("                      ______________\n");
    cat("                     |    ddsPLS    |\n");
    cat("=====================----------------=====================\n");
  }
  while (test){
    if (verbose) {
      cat(paste("Should we build component " ,h+1 , " ? Bootstrap pending...\n",sep=""))
    }
    NCORES_w <- min(NCORES,n_B)
    n_B_i <- ceiling(n_B/NCORES)
    `%my_do%` <- ifelse(NCORES_w!=1,{
      out<-`%dopar%`;cl <- makeCluster(NCORES_w)
      registerDoParallel(cl);out},{out <- `%do%`;out})
    res <- foreach(i_B=1:NCORES_w,.packages = "testEigen",
                   .combine='c',.multicombine=TRUE) %my_do% {
                     bootstrapWrap(U_out,V0,X_init,Y_init,lambdas,lambda_prev,
                                   R=h+1,n_B_i,doBoot=TRUE,n,p,q,N_lambdas)
                   }
    if(NCORES_w!=1) stopCluster(cl)
    Results$R2[[h+1]] <- do.call(rbind,res[which(names(res)=="R2")])
    Results$R2h[[h+1]] <- do.call(rbind,res[which(names(res)=="R2h")])
    Results$Q2[[h+1]] <- do.call(rbind,res[which(names(res)=="Q2")])
    Results$Q2h[[h+1]] <- do.call(rbind,res[which(names(res)=="Q2h")])
    Results$PropQ2hPos[[h+1]] <- apply(Results$Q2h[[h+1]],2,function(v){sum(v>=0)/length(v)})
    Results$R2mean[[h+1]] <- colMeans(Results$R2[[h+1]])
    Results$R2hmean[[h+1]] <- colMeans(Results$R2h[[h+1]])
    Results$Q2mean[[h+1]] <- colMeans(Results$Q2[[h+1]])
    Results$Q2hmean[[h+1]] <- colMeans(Results$Q2h[[h+1]])
    Results$R2mean_diff_Q2mean[[h+1]] <- Results$R2mean[[h+1]]-Results$Q2mean[[h+1]]
    TEST <- (Results$Q2hmean[[h+1]]>lowExplainedVariance)*
      (Results$Q2mean[[h+1]]>Q2_previous)*
      (Results$PropQ2hPos[[h+1]]>minBootProp)==1
    nb_ValsOk = sum(TEST)
    test_lambdas[[h+1]] <- TEST
    # resOUT <- NULL
    test_t2 <- F
    if (nb_ValsOk>0){
      bestVal = min(Results$R2mean_diff_Q2mean[[h+1]][TEST])
      bestID = which(Results$R2mean_diff_Q2mean[[h+1]]==bestVal)[1]
      lambda_prev[h+1] = lambdas[bestID]
      resMozna <- modelddsPLSCpp_Rcpp(U_out,V0,X_init,Y_init,lambda_prev,R=h+1,n,p,q)
      test_t2 <- sum((resMozna$t[,h+1])^2)>errorMin
      if(test_t2){
        resMozna -> resOUT
        resMozna <- NULL
        h <- h + 1
        Q2Best[h] = Results$Q2mean[[h]][bestID]
        Q2hBest[h] = Results$Q2hmean[[h]][bestID]
        R2Best[h] = Results$R2mean[[h]][bestID]
        R2hBest[h] = Results$R2hmean[[h]][bestID]
        Q2_previous = Q2Best[h]
        U_out[,1:h] = resOUT$U[,1:h]
        V0[,1:h] = resOUT$V[,1:h]
        B_out <- resOUT$B
        for (i in 1:p){
          if(sdX[i]>errorMin){
            B_out[i,] <- B_out[i,]/sdX[i]
          }
        }
        for (j in 1:q){
          B_out[,j] <- B_out[,j]*sdY[j]
        }
        B_out -> resOUT$B
        out0 <- list(model=list(muX=muX,muY=muY,B=B_previous));class(out0)="ddsPLS"
        out1 <- list(model=list(muX=muX,muY=muY,B=B_out));class(out1)="ddsPLS"
        Y_est_0 <- predict(out0,X);ind_1 <- (colSums(abs(B_out))>1e-9)*1
        Y_est_1 <- predict(out1,X);ind_0 <- (colSums(abs(B_previous))>1e-9)*1
        cor2_0 <- unlist(lapply(1:q,function(j){if(ind_0[j]==1){
          1-sum((Y_est_0[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
          }else{0}
          }))
        cor2_1 <- unlist(lapply(1:q,function(j){if(ind_1[j]==1){
          1-sum((Y_est_1[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
          }else{0}}))
        B_previous <- B_out
        varExplained <- c(varExplained,mean(cor2_1-cor2_0)*100)
        varExplainedTot <- c(varExplainedTot,mean(cor2_1)*100)
        varExplained_y <- rbind(varExplained_y,(cor2_1-cor2_0)*100)
        varExplainedTot_y <- rbind(varExplainedTot_y,(cor2_1)*100)
        if (verbose) {
          ress <- data.frame(
            list(
              "  "="   ",
              "lambda"=round(lambdas[bestID],2),
              "R2"=round(Results$R2mean[[h]][bestID],2),
              "R2h"=round(Results$R2hmean[[h]][bestID],2),
              "Q2"=round(Results$Q2mean[[h]][bestID],2),
              "Q2h"=round(Results$Q2hmean[[h]][bestID],2),
              "VarExpl"=paste(round(varExplained[h]),"%",sep=""),
              "VarExpl Tot"=paste(round(varExplainedTot[h]),"%",sep="")
            )
          )
          rownames(ress) <- ""
          colnames(ress)[1] <- "   "
          print(ress)
          cat(paste("                                     ...component ",h," built!","\n",sep=""))
        }
      }
    }
    if (!test_t2 | nb_ValsOk<=0){
      if(verbose) cat(paste("                                 ...component ",h+1," not built!","\n",sep=""))
      test = F;
      if(h==0){
        if(verbose){
          if(sum(Results$Q2hmean[[h+1]]>lowExplainedVariance)==0){
            cat("             ...no Q2r large enough for tested lambda.\n")
          }
        }
      }
    }
  }
  if(verbose) {
    cat("=====================                =====================\n");
    cat("                     ================\n");
  }
  lambda_sol=R2Sol=R2hSol=Q2Sol=Q2hSol <- rep(0,h)
  for (r in 0:h) {
    lambda_sol[r] = lambda_prev[r];
    R2Sol[r] = R2Best[r];
    R2hSol[r] = R2hBest[r];
    Q2Sol[r] = Q2Best[r];
    Q2hSol[r] = Q2hBest[r];
  }
  out <- list()
  if(h>0){
    out$model <- resOUT
    out$model$muY <- muY
    out$model$muX <- muX
    out$model$sdY <- sdY
    out$model$sdX <- sdX
    Results$lambdas <- lambdas
    out$results <- Results
    out$varExplained <- list()
    out$varExplained$Comp <- varExplained
    out$varExplained$Cumu <- varExplainedTot
  }else{
    out$model = NULL
    out$results <- NULL
  }
  out$R = h
  out$lambda = lambda_sol
  out$lambda_optim <- test_lambdas
  out$Q2 = Q2Sol
  out$Q2h = Q2hSol
  out$R2 = R2Sol
  out$R2h = R2hSol
  out$lowExplainedVariance=lowExplainedVariance
  class(out) <- "ddsPLS"
  if(h>0){
    out$Y_est <- predict(out,X)
    out$Y_obs <- Y
    colnames(out$Y_est) = colnames(Y)
    rownames(out$model$V) = colnames(Y)
    rownames(out$model$U) = colnames(X)
    rownames(out$model$U_star) = colnames(X)
    rownames(out$model$B) = colnames(X)
    colnames(out$model$B) = colnames(Y)
    out$varExplained$PerY <- varExplainedTot_y[h,]
    out$varExplained$PerYPerComp <- list()
    out$varExplained$PerYPerComp$Comp <- varExplained_y
    out$varExplained$PerYPerComp$Cumu <- varExplainedTot_y
    selX <- (rowSums(abs(out$model$B))>1e-9)*1
    selY <- (colSums(abs(out$model$B))>1e-9)*1
    out$Selection <- list(X=which(selX==1),Y=which(selY==1))
  }
  if (verbose & h>0) {
    plot(out)
  }
  return(out)
}
