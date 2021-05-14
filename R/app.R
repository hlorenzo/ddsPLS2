#' Applet to start ddsPLS
#'
#' @param ...
#'
#' @return
#' @export
#'
ddsPLS2_App <- function(...) {
  vizu <- c("predict","Q2","criterion", "Q2r", "R2r", "R2", "weightX", "weightY","loadingX","loadingY")
  pospos <- c("topleft","topright","bottomright","bottomleft",
              "center","top","bottom",
              "left","right")
  lamCh <- c("lambda.min","lambda.1se")
  get_datas <- function(n=200,p2=100,sd2=3){
    phi <- matrix(rnorm(2*n),nrow = n)
    y <- phi[,1,drop=FALSE] + rnorm(n,sd = 0.3)
    p1_1 <- 50
    p1_2 <- 100
    p1_3 <- 50
    p2_1 <- 100
    x1 <- cbind(matrix(rep(phi[,1,drop=FALSE],p1_1),byrow = F,nrow = n) + rnorm(n*p1_1,sd = 0.4),
                matrix(rep(phi[,1,drop=FALSE]+phi[,2,drop=FALSE],p1_2),byrow = F,nrow = n) + rnorm(n*p1_2,sd = 0.4),
                matrix(rep(phi[,2,drop=FALSE],p1_3),byrow = F,nrow = n) + rnorm(n*p1_3,sd = 0.4))
    x2<- cbind(matrix(rep(phi[,1,drop=FALSE],p2_1),byrow = F,nrow = n) + rnorm(n*p2_1,sd = sd2),
               matrix(rnorm(n*p2,sd=sd2),byrow = F,nrow = n))
    list(x1,x2,y)
  }
  get_datas_NoX1 <- function(n=200,p2=100,sd2=3){
    phi <- matrix(rnorm(2*n),nrow = n)
    y <- phi[,1,drop=FALSE] + rnorm(n,sd = 0.3)
    p1_1 <- 50
    p1_2 <- 100
    p1_3 <- 50
    p2_1 <- 100
    x1 <- cbind(matrix(rep(phi[,1,drop=FALSE],p1_1),byrow = F,nrow = n) + rnorm(n*p1_1,sd = 0.4),
                matrix(rep(phi[,1,drop=FALSE]+phi[,2,drop=FALSE],p1_2),byrow = F,nrow = n) + rnorm(n*p1_2,sd = 0.4),
                matrix(rep(phi[,2,drop=FALSE],p1_3),byrow = F,nrow = n) + rnorm(n*p1_3,sd = 0.4))
    x2<- cbind(matrix(rep(phi[,1,drop=FALSE],p2_1),byrow = F,nrow = n) + rnorm(n*p2_1,sd = sd2),
               matrix(rnorm(n*p2,sd=sd2),byrow = F,nrow = n))
    list(x1[,-c(1:p1_1)],x2,y)
  }
  cols_gps <- RColorBrewer::brewer.pal(5,"Set1")
  ui <- fluidPage(
    #=======================================
    navbarPage("ddsPLS (data-driven Sparse PLS)",
               tabPanel("Uploading Files",
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("fileX", "Choose CSV Files for X",multiple = TRUE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                            fileInput("fileY", "Choose CSV File for Y",multiple = TRUE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                            tags$hr(),
                            checkboxInput("header", "Header", TRUE),
                            radioButtons("sep", "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ","),
                            radioButtons("dec", "Decimal",choices = c(Point = ".",Comma = ","),selected = "."),
                            radioButtons("quote", "Quote",choices = c(None = "","Double Quote" = '"',"Single Quote" = "'"),selected = '"'),
                            tags$hr(),
                            actionButton("startimport","Upload",icon=icon("play")),
                            numericInput("sizeplot2", "Size of plot (pixels)",600,min=200,max=3000,step = 100),
                            tableOutput("contents")
                          ),
                          mainPanel(
                            plotOutput('plot1R')
                          )
                        )),
               tabPanel("Building model",
                        sidebarLayout(
                          sidebarPanel(
                            numericInput('n_B', 'Number of Bootstrap samples',50,min=50,max=1000,step = 50),
                            actionButton("run","Run analysis",icon=icon("play"))
                          ),
                          mainPanel(
                            verbatimTextOutput("summary")
                          )
                        )
               ),
               tabPanel("Vizualize results",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput('plo', 'Type of vizualization', vizu),
                            selectInput('pos', 'Legend position', pospos),
                            numericInput('sizeplot', 'Size of plot (pixels)',600,min=200,max=3000,step = 100)
                          ),
                          mainPanel(
                            plotOutput('plot2')
                          )
                        )
               ),
               tabPanel("ddsPLS VS PLS",
                        sidebarLayout(
                          sidebarPanel(
                            numericInput('ncomp', 'Number of components for PLS model', 1,min=1,max=10,step = 1)
                          ),
                          mainPanel(
                            plotOutput('plot3')
                          )
                        )
               )
    )
    #=======================================
  )
  server <- function(input, output, session) {
    #=======================================
    #=======================================
    datasR <- eventReactive(input$startimport, {
      K <- nrow(input$fileX)
      tryCatch(
        {
          dfX_list <- lapply(input$fileX$datapath,function(fifi){
            read.csv(fifi,header = input$header,sep = input$sep,dec = input$dec,quote = input$quote)
          })
          dfY <- read.csv(input$fileY$datapath,header = input$header,sep = input$sep,dec = input$dec,quote = input$quote)
        },
        error = function(e) {
          stop(safeError(e))
        }
      )
      list(isSimu=F,Xs = dfX_list,Y=dfY,
           ps=unlist(lapply(dfX_list,ncol)),
           colsReal=unlist(lapply(1:length(ps),function(k){rep(k,ps[k])})))
    })

    output$contents <- eventReactive(input$startimport, {
      dada <- datasR()
      dfX_list <- dada$Xs
      dfY <- dada$Y
      K <- length(dfX_list)
      out <- matrix(NA,K+1,3)
      colnames(out) <- c("Block","Number of rows","Number of columns")
      for(k in 1:K){
        out[k,2] <- nrow(dfX_list[[k]])
        out[k,3] <- ncol(dfX_list[[k]])
      }
      out[K+1,2] <- nrow(dfY)
      out[K+1,3] <- ncol(dfY)
      out[,1] <- c(paste("X",1:K),"Y")
      return(data.frame(out))
    })

    sizeplot2 <- reactive({
      input$sizeplot2
    })

    output$plot1R <- renderPlot({
      K <- length(datasR()$Xs)
      par(mfrow=c(1,K),mar=c(2,3,3,1))
      for(k in 1:K){
        matplot((cor(as.matrix(datasR()$Xs[[k]]), as.matrix(datasR()$Y))),type="l",
                main=paste("Correlation y/Block",k),ylim=c(-1,1),ylab="",lty=1)
        abline(h=c(0,1,-1))
      }
    },width = sizeplot2)

    model <- eventReactive(input$run, {
      req(input$fileX,input$fileY)
      x <- as.matrix(do.call(cbind,datasR()$Xs))
      y <- as.matrix(datasR()$Y)
      mo <- ddsPLS(x,y,
                   verbose=F,doBoot = T,
                   lambdas = NULL,n_B = input$n_B)
      return(mo)
    })

    output$summary <- renderPrint({
      summary(model())
    })

    sizeplot <- reactive({
      input$sizeplot
    })

    output$plot2 <- renderPlot({
      mo <- model()
      noModel <- mo$R==0
      if(!noModel){
        colo <- datasR()$colsReal
        plot(mo,type = input$plo,legend.position =input$pos,col=colo)
      }else{
        plot(0,0,col="white",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
        text(x = 0,y=0,"Nothing to be plotted because model is empty.")
      }
    },height = sizeplot)

    bPLS <- reactive({
      req(input$fileX,input$fileY)
      x <- as.matrix(do.call(cbind,datasR()$Xs))
      y <- as.matrix(datasR()$Y)
      df <- data.frame(x=x,y=y)
      moPLS <- plsr(y~.,data=df,ncomp=input$ncomp)
      out <- moPLS$coefficients[,,input$ncomp]
      if(!is.null(out)){
        if(is.vector(out)){
          out <- matrix(out,ncol=1)
        }
      }
      return(out)
    })

    output$plot3 <- renderPlot({
      req(input$fileX,input$fileY)
      dada <- datasR()
      colo <- dada$colsReal
      ps <- dada$ps
      K <- length(ps)
      q <- ncol(dada$Y)
      # }
      p_sum <- c(0,cumsum(ps))
      layout(matrix(1:(2*K),nrow = 2,byrow = F))
      for(k in 1:K){
        id_var <- c(p_sum[k]+1):p_sum[k+1]
        bPLS <- bPLS()[id_var,]
        bddsPLS <- model()$model$B[id_var,]
        matplot(bPLS,type="l",col="gray",ylim=c(-1,1)*max(abs(bPLS)),lty=1,
                main=paste("PLS regression coefficients, block X",k),xlab="Index",ylab="")
        abline(h=0)
        matplot(bPLS,type="p",add=T,xlab="",ylab="",col=1:q,pch=1:q)
        matplot(bddsPLS,type="l",col="gray",lty=1,ylim=c(-1,1)*max(abs(bddsPLS)),
                main=paste("ddsPLS regression coefficients, block X",k),xlab="Index",ylab="")
        abline(h=0)
        matplot(bddsPLS,type="p",add=T,xlab="",ylab="",pch=1:q,col=1:q)
      }
    },height = 1000)
    #=======================================
    #=======================================
  }
  shinyApp(ui, server, ...)
}