#' @title Linear Adsorption Kinetic Model Summary
#' @description Summarized results of parameter and error values collected from linear adsorption kinetic models, namely: Elovich, Fractional Power, Pseudo-First-Order, and Pseudo-Second-Order.
#' @param t the numerical value for contact time. This parameter should not be equal to zero to prevent infinite value. Any row(s) that contain(s) value of t equal to zero will be automatically removed to proceed with the calculation.
#' @param qt the numerical value for the amount adsorbed at time t. This parameter should not be equal to qe or zero as it will cause an infinite value. Any row(s) that contain(s) value of qt equal to qe or zero will be automatically removed to proceed with the calculation.
#' @param qe the numerical value for the amount adsorbed at equilibrium
#' @param sort.by the name of the statistical error parameter in which the models are sorted in either increasing or decreasing order. The only accepted arguments are "RMSE" for Relative Mean Square Error, 'MAE' for Mean Absolute Error, 'MSE' for Mean Squared Error, 'RAE' for Relative Absolute Error, 'AIC' for Akaike Information Criterion, 'BIC' for Bayesian Information Criterion, 'R2' for Coefficient of Determination, and 'SE' for Standard Error Estimate. This argument is case-sensitive, and failure to input the correct value will yield a summary of models in alphabetical order.
#' @import stats
#' @import Metrics
#' @import utils
#' @return the summarized error and parameter values from adsorption kinetic models.
#' @examples t  <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples qe <- 4.8
#' @examples aksummary.l(t,qt,qe,"R2")
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @export

aksummary.l <- function(t,qt,qe,sort.by){
  errors <- new.env()
  parameters <- new.env()
  summary.data <- new.env()
  if(missing(qe)){
    stop("qe is required to run the function. The qe value can be solve using isotherm models.
         See pakcage PUPAIM to solve adsorption isotherm models")
  }else{qe <- qe}
  t <-t ;qt <-qt
  if(missing(sort.by)){
    s<-"sort.by"
  }else if(is.null(sort.by)){
    s<-"sort.by"
  }else if(isFALSE(sort.by)){
    s<-"sort.by"
  }else{
    s<-sort.by
  }
  PFO.sum.l<- function(t,qt,qe){
    qe <- qe         ;x   <- t
    y  <- log(qe-qt) ;int <- log(qe)
    dat1  <- data.frame(x,y)
    dat   <- dat1[complete.cases(dat1),]
    n1    <- nrow(dat1)
    n.dat <- nrow(dat)
    y <- dat$y       ;x <- dat$x
    fit37 <- lm(I(y -int) ~ 0+x)
    lin.val <- coef(fit37)
    slp <- as.numeric(lin.val[1])
    pred.val <- (int+(slp*x))
    summary.data$PFO.l<-summary(fit37)
    parameters$pfo.l.k1 <- -slp
    errors$rmse.pfo.l <- rmse(y,pred.val)
    errors$mae.pfo.l  <- mae(y,pred.val)
    errors$mse.pfo.l  <- mse(y,pred.val)
    errors$rae.pfo.l  <- rae(y,pred.val)
    errors$PAIC.pfo.l <- AIC(fit37)
    errors$PBIC.pfo.l <- BIC(fit37)
    errors$SE.pfo.l   <- sqrt((sum((y-pred.val)^2))/(n.dat-2))
    errors$R2.pfo.l   <- summary(fit37)$r.squared
  }
  PSO.sum.l<- function(t,qt,qe){
    x  <- t          ;y <- t/qt
    qe <- qe
    dat1  <- data.frame(x,y)
    dat   <- dat1[complete.cases(dat1),]
    n1    <- nrow(dat1)
    n.dat <- nrow(dat)
    y <- dat$y       ;x <- dat$x
    slp <- (1/qe)
    fit38 <- lm(y - slp*x ~ 1)
    lin.val <- coef(fit38)
    int<-as.numeric(lin.val[1])
    pred.val <- ((slp*x)+int)
    RSS <- sum((y-pred.val)^2)
    TSS <- sum((y-mean(y))^2)
    summary.data$PSO.l<-summary(fit38)
    parameters$pso.l.k2<- (1/((qe^2)*int))
    errors$rmse.pso.l <- rmse(y,pred.val)
    errors$mae.pso.l  <- mae(y,pred.val)
    errors$mse.pso.l  <- mse(y,pred.val)
    errors$rae.pso.l  <- rae(y,pred.val)
    errors$PAIC.pso.l <- AIC(fit38)
    errors$PBIC.pso.l <- BIC(fit38)
    errors$SE.pso.l   <- sqrt((sum((y-pred.val)^2))/(n.dat-2))
    errors$R2.pso.l   <- 1-(RSS/TSS)
  }
  Elovich.sum.l<- function(t,qt){
    x <- log(t)  ;y <- qt
    dat1 <- data.frame(x,y)
    dat  <- subset(dat1, x!="-Inf")
    n1      <- nrow(na.omit(dat1))
    n.dat   <- nrow(na.omit(dat))
    y <- dat$y   ;x <- dat$x
    fit39 <- lm(y ~ x)
    lin.val <- coef(fit39)
    int <-as.numeric(lin.val[1])
    slp <-as.numeric(lin.val[2])
    pred.val <- (int+(slp*x))
    summary.data$Elovich.l <-summary(fit39)
    parameters$elovich.l.b <- 1/slp
    parameters$elovich.l.a <- exp(int*parameters$elovich.l.b) + parameters$elovich.l.b
    errors$rmse.elovich.l <- rmse(y,pred.val)
    errors$mae.elovich.l  <- mae(y,pred.val)
    errors$mse.elovich.l  <- mse(y,pred.val)
    errors$rae.elovich.l  <- rae(y,pred.val)
    errors$PAIC.elovich.l <- AIC(fit39)
    errors$PBIC.elovich.l <- BIC(fit39)
    errors$SE.elovich.l   <- sqrt((sum((y-pred.val)^2))/(n.dat-2))
    errors$R2.elovich.l   <- cor(x,y)^2
  }
  FractionalPower.sum.l<- function(t,qt,qe){
    qe <- qe       ;x <- log(t)
    y  <- log(qt/qe)
    dat1  <- data.frame(x,y)
    dat   <- subset(dat1, x!="-Inf")
    dat   <- subset(dat,  y!="-Inf")
    n1    <- nrow(na.omit(dat1))
    n.dat <- nrow(na.omit(dat))
    y  <- dat$y    ;x <- dat$x
    fit40 <- lm(y ~ x)
    lin.val <- coef(fit40)
    int<-as.numeric(lin.val[1])
    slp <-as.numeric(lin.val[2])
    pred.val <- (int+(slp*x))
    summary.data$FractionalPower.l<-summary(fit40)
    parameters$FractionalPower.l.b <- slp
    parameters$FractionalPower.l.a <- exp(int)
    errors$rmse.FractionalPower.l <- rmse(y,pred.val)
    errors$mae.FractionalPower.l  <- mae(y,pred.val)
    errors$mse.FractionalPower.l  <- mse(y,pred.val)
    errors$rae.FractionalPower.l  <- rae(y,pred.val)
    errors$PAIC.FractionalPower.l <- AIC(fit40)
    errors$PBIC.FractionalPower.l <- BIC(fit40)
    errors$SE.FractionalPower.l   <- sqrt((sum((y-pred.val)^2))/(n.dat-2))
    errors$R2.FractionalPower.l   <- cor(x,y)^2
  }

  Elovich.sum.l(t,qt)
  FractionalPower.sum.l(t,qt,qe)
  PFO.sum.l(t,qt,qe)
  PSO.sum.l(t,qt,qe)
  Elovich.alpha <- parameters[["elovich.l.a"]]
  Elovich.beta  <- parameters[["elovich.l.b"]]
  FractionalPower.alpha <- parameters[["FractionalPower.l.a"]]
  FractionalPower.beta  <- parameters[["FractionalPower.l.b"]]
  PFO.k1 <- parameters[["pfo.l.k1"]]
  PSO.k2 <- parameters[["pso.l.k2"]]
  param1.name <- c("alpha","alpha","k.PFO","k.PSO")
  param2.name <- c("beta","beta","-","-")
  param1.val  <- as.numeric(c(Elovich.alpha,FractionalPower.alpha,PFO.k1,PSO.k2))
  param2.val  <- as.numeric(c(Elovich.beta,FractionalPower.beta," "," "))
  param1.val  <- round(param1.val,digits=6)
  param2.val  <- round(param2.val,digits=6)
  param2.val[is.na(param2.val)]<-"-"
  Models <- c("Elovich","FractionalPower","PFO","PSO")
  rmse.val <- as.numeric(c(errors[["rmse.elovich.l"]],errors[["rmse.FractionalPower.l"]],errors[["rmse.pfo.l"]],errors[["rmse.pso.l"]]))
  mae.val  <- as.numeric(c(errors[["mae.elovich.l"]],errors[["mae.FractionalPower.l"]],errors[["mae.pfo.l"]],errors[["mae.pso.l"]]))
  mse.val  <- as.numeric(c(errors[["mse.elovich.l"]],errors[["mse.FractionalPower.l"]],errors[["mse.pfo.l"]],errors[["mse.pso.l"]]))
  rae.val  <- as.numeric(c(errors[["rae.elovich.l"]],errors[["rae.FractionalPower.l"]],errors[["rae.pfo.l"]],errors[["rae.pso.l"]]))
  PAIC.val <- as.numeric(c(errors[["PAIC.elovich.l"]],errors[["PAIC.FractionalPower.l"]],errors[["PAIC.pfo.l"]],errors[["PAIC.pso.l"]]))
  PBIC.val <- as.numeric(c(errors[["PBIC.elovich.l"]],errors[["PBIC.FractionalPower.l"]],errors[["PBIC.pfo.l"]],errors[["PBIC.pso.l"]]))
  R2.val   <- as.numeric(c(errors[["R2.elovich.l"]],errors[["R2.FractionalPower.l"]],errors[["R2.pfo.l"]],errors[["R2.pso.l"]]))
  SE.val   <- as.numeric(c(errors[["SE.elovich.l"]],errors[["SE.FractionalPower.l"]],errors[["SE.pfo.l"]],errors[["SE.pso.l"]]))
  rmse.val <- round(rmse.val,digits=6)
  mae.val  <- round(mae.val,digits=6)
  mse.val  <- round(mse.val,digits=6)
  rae.val  <- round(rae.val,digits=6)
  PAIC.val <- round(PAIC.val,digits=4)
  PBIC.val <- round(PBIC.val,digits=4)
  R2.val   <- round(R2.val,digits=8)
  SE.val   <- round(SE.val,digits=6)
  Summary.Table <- data.frame(Models,rmse.val,mae.val,mse.val,rae.val,PAIC.val,PBIC.val,R2.val,SE.val,param1.name,param1.val,param2.name,param2.val)
  colnames(Summary.Table) <- c("Models","RMSE","MAE","MSE","RAE","AIC","BIC","R2","SE","Parameter1","     ","Parameter2","      ")
  if(s=="RMSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RMSE),]
    top.model <- as.character(Summary.Table$Models[which.min(rmse.val)])
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Root Mean Square Error (RMSE)")
  }else if(s=="MAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MAE),]
    top.model <- as.character(Summary.Table$Models[which.min(mae.val)])
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Mean Absolute Error (MAE)")
  }else if(s=="MSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MSE),]
    top.model <- as.character(Summary.Table$Models[which.min(mse.val)])
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Mean Squared Error (MSE)")
  }else if(s=="RAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RAE),]
    top.model <- as.character(Summary.Table$Models[which.min(rae.val)])
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Relative Absolute Error (RAE)")
  }else if(s=="AIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$AIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PAIC.val)])
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Akaike Information Criterion (AIC)")
  }else if(s=="BIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$BIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PBIC.val)])
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Bayesian Information Criterion (BIC)")
  }else if(s=="R2"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$R2,decreasing=TRUE),]
    top.model <- as.character(Summary.Table$Models[which.max(R2.val)])
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Coefficient of Determination (R2)")
  }else if(s=="SE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$SE),]
    top.model <- as.character(Summary.Table$Models[which.min(SE.val)])
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Standard Error (SE)")
  }else{
    Sort.Summary.Table <- Summary.Table
    message("Summary of Linear Adsorption Kinetic Models: List of Error Values and Parameters (Unsorted)")
    top.model <- as.character("Unsorted")
  }
  Param.Table <- Sort.Summary.Table[c(-2,-3,-4,-5,-6,-7,-8,-9)]
  Sort.Summary.Table <-Sort.Summary.Table[c(-10,-11,-12,-13)]
  print(Sort.Summary.Table, right=FALSE,row.names = FALSE)
  print(Param.Table, right=FALSE,row.names = FALSE )
  if(as.character(top.model)!= "Unsorted"){
    message("Summary of the most fitted model (",top.model,")",sep="")
    top.model<-capture.output(message(top.model,".l",sep=""),type="message")
    print(summary.data[[top.model]])
  }
}
