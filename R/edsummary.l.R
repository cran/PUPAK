#' @title Linear External Diffusion Models Summary
#' @description Summarized results of parameter and error values collected from external diffusion models, namely: Boyd External Diffusion, Furusawa and Smith, and Matthews and Weber.
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param Ct the numerical value for the concentration of the adsorbent at time t
#' @param qinf the numerical value for the amount adsorbed at infinite time
#' @param Co the numerical value for the initial concentration of the adsorbent
#' @param m the numerical value for mass of adsorbent
#' @param V the numerical value for volume of solution
#' @param b the numerical value or the Langmuir isotherm constant
#' @param sort.by the name of the statistical error parameter in which the models are sorted in either increasing or decreasing order. The only accepted arguments are "RMSE" for Relative Mean Square Error, 'MAE' for Mean Absolute Error, 'MSE' for Mean Squared Error, 'RAE' for Relative Absolute Error, 'AIC' for Akaike Information Criterion, 'BIC' for Bayesian Information Criterion, 'R2' for Coefficient of Determination, and 'SE' for Standard Error Estimate. This argument is case-sensitive, and failure to input the correct value will yield a summary of models in alphabetical order.
#' @import stats
#' @import Metrics
#' @import utils
#' @return the summarized error and parameter values from external diffusion models.
#' @examples t  <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples Ct <- c(10.000,8.141,8.056,7.949,7.863,7.799,7.778,7.756,7.692)
#' @examples edsummary.l(t,qt,Ct,qinf=4.8,Co=10,m=0.05,V=0.1,b=1.3,"SE")
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @export
edsummary.l <- function(t,qt,Ct,qinf,Co,m,V,b,sort.by){
  errors <- new.env()
  parameters <- new.env()
  summary.data <- new.env()
  if(missing(qinf)){
    stop("qinf is required to run the function")
  }else if(missing(Co)){
    stop("Co is required to run the function")
  }else if(missing(m)){
    stop("m is required to run the function")
  }else if(missing(V)){
    stop("V is required to run the function")
  }else if(missing(b)){
    stop("b is required to run the function. the parameter b is equal to the Langmuir isotherm parameter.
         See pakcage PUPAIM to solve adsorption isotherm models")
  }else{}
  t<-t     ; b <- b
  qt<-qt   ; qinf <-qinf
  Ct <- Ct ; Co <- Co
  m <- m   ; V <- V
  if(missing(sort.by)){
    s<-"sort.by"
  }else if(is.null(sort.by)){
    s<-"sort.by"
  }else if(isFALSE(sort.by)){
    s<-"sort.by"
  }else{
    s<-sort.by
  }
  BED.sum.l<- function(t,qt,qinf){
    qinf <- qinf
    x    <- t
    y    <- log(1-(qt/qinf))
    dat1 <- data.frame(x,y)
    dat  <- subset(dat1,  y!="-Inf")
    dat  <- subset(dat,  y!="NaN")
    n1   <- nrow(dat1)
    n.dat<- nrow(dat)
    y <- dat$y   ;x <- dat$x
    fit41 <- lm(y ~ x)
    lin.val <- coef(fit41)
    int <-as.numeric(lin.val[1])
    slp <-as.numeric(lin.val[2])
    pred.val <- (int+(slp*x))
    summary.data$Boyd.External.Diffusion.l<-summary(fit41)
    parameters$BED.l.R <- -slp
    parameters$BED.l.A <- int
    errors$rmse.BED.l <- rmse(y,pred.val)
    errors$mae.BED.l  <- mae(y,pred.val)
    errors$mse.BED.l  <- mse(y,pred.val)
    errors$rae.BED.l  <- rae(y,pred.val)
    errors$PAIC.BED.l <- AIC(fit41)
    errors$PBIC.BED.l <- BIC(fit41)
    errors$SE.BED.l   <- sqrt((sum((y-pred.val)^2))/(n.dat-2))
    errors$R2.BED.l   <- cor(x,y)^2
  }
  FAS.sum.l<- function(t,Ct,m,V,b,Co){
    m <- m ;V  <- V  ;b <- b
    x <- t ;Co <- Co
    y <- (1/(1+(1/((m/V)*b))))*log((Ct/Co)-((1/((m/V)*b))*(1-(Ct/Co))))
    dat1   <- data.frame(x,y)
    dat    <- subset(dat1,  y!="-Inf")
    dat    <- subset(dat,  y!="NaN")
    n1     <- nrow(dat1)
    n.dat  <- nrow(dat)
    y   <- dat$y   ;x   <- dat$x
    fit42 <- lm(y ~ x)
    lin.val <- coef(fit42)
    int<-as.numeric(lin.val[1])
    slp <-as.numeric(lin.val[2])
    pred.val <- (int+(slp*x))
    summary.data$Furusawa.and.Smith.l<-summary(fit42)
    parameters$FAS.l.kFASS <- -slp
    errors$rmse.FAS.l <- rmse(y,pred.val)
    errors$mae.FAS.l  <- mae(y,pred.val)
    errors$mse.FAS.l  <- mse(y,pred.val)
    errors$rae.FAS.l  <- rae(y,pred.val)
    errors$PAIC.FAS.l <- AIC(fit42)
    errors$PBIC.FAS.l <- BIC(fit42)
    errors$SE.FAS.l   <- sqrt((sum((y-pred.val)^2))/(n.dat-2))
    errors$R2.FAS.l   <- cor(x,y)^2
  }
  MAW.sum.l<- function(t,Ct,Co){
    Co <- Co
    x  <- t
    y  <- log(Ct/Co)
    dat1 <- data.frame(x,y)
    dat  <- subset(dat1,  y!="-Inf")
    dat  <- subset(dat,  y!="NaN")
    n1    <- nrow(dat1)
    n.dat   <- nrow(dat)
    y <- dat$y  ;x <- dat$x
    fit43 <- lm(y ~ x)
    lin.val <- coef(fit43)
    int<-as.numeric(lin.val[1])
    slp <-as.numeric(lin.val[2])
    pred.val <- (int+(slp*x))
    summary.data$Matthews.and.Weber.l<-summary(fit43)
    parameters$MAW.l.kMAWS <- -slp
    errors$rmse.MAW.l <- rmse(y,pred.val)
    errors$mae.MAW.l  <- mae(y,pred.val)
    errors$mse.MAW.l  <- mse(y,pred.val)
    errors$rae.MAW.l  <- rae(y,pred.val)
    errors$PAIC.MAW.l <- AIC(fit43)
    errors$PBIC.MAW.l <- BIC(fit43)
    errors$SE.MAW.l   <- sqrt((sum((y-pred.val)^2))/(n.dat-2))
    errors$R2.MAW.l   <- cor(x,y)^2
  }

  BED.sum.l(t,qt,qinf)
  FAS.sum.l(t,Ct,m,V,b,Co)
  MAW.sum.l(t,Ct,Co)
  BED.R <- parameters[["BED.l.R"]]
  BED.A <- parameters[["BED.l.A"]]
  FAS.kFASS <- parameters[["FAS.l.kFASS"]]
  MAW.kMAWS <- parameters[["MAW.l.kMAWS"]]
  param1.name <- c("R","kFASS","kMAWS")
  param2.name <- c("A","-","-")
  param1.val <- as.numeric(c(BED.R,FAS.kFASS,MAW.kMAWS))
  param2.val <- as.numeric(c(BED.A," "," "))
  param1.val <- round(param1.val,digits=6)
  param2.val <- round(param2.val,digits=6)
  param2.val[is.na(param2.val)]<-"-"
  Models <- c("Boyd.External.Diffusion","Furusawa.and.Smith","Matthews.and.Weber")
  rmse.val <- as.numeric(c(errors[["rmse.BED.l"]],errors[["rmse.FAS.l"]],errors[["rmse.MAW.l"]]))
  mae.val  <- as.numeric(c(errors[["mae.BED.l"]],errors[["mae.FAS.l"]],errors[["mae.MAW.l"]]))
  mse.val  <- as.numeric(c(errors[["mse.BED.l"]],errors[["mse.FAS.l"]],errors[["mse.MAW.l"]]))
  rae.val  <- as.numeric(c(errors[["rae.BED.l"]],errors[["rae.FAS.l"]],errors[["rae.MAW.l"]]))
  PAIC.val <- as.numeric(c(errors[["PAIC.BED.l"]],errors[["PBIC.FAS.l"]],errors[["PBIC.MAW.l"]]))
  PBIC.val <- as.numeric(c(errors[["PAIC.BED.l"]],errors[["PAIC.FAS.l"]],errors[["PAIC.MAW.l"]]))
  R2.val   <- as.numeric(c(errors[["R2.BED.l"]],errors[["R2.FAS.l"]],errors[["R2.MAW.l"]]))
  SE.val   <- as.numeric(c(errors[["SE.BED.l"]],errors[["SE.FAS.l"]],errors[["SE.MAW.l"]]))
  rmse.val <- round(rmse.val,digits=6)
  mae.val  <- round(mae.val,digits=6)
  mse.val  <- round(mse.val,digits=6)
  rae.val  <- round(rae.val,digits=6)
  PAIC.val <- round(PAIC.val,digits=4)
  PBIC.val <- round(PBIC.val,digits=4)
  R2.val   <- round(R2.val,digits=4)
  SE.val   <- round(SE.val,digits=6)
  Summary.Table <- data.frame(Models,rmse.val,mae.val,mse.val,rae.val,PAIC.val,PBIC.val,R2.val,SE.val,param1.name,param1.val,param2.name,param2.val)
  colnames(Summary.Table) <- c("Models","RMSE","MAE","MSE","RAE","AIC","BIC","R2","SE","Parameter1","     ","Parameter2","      ")
  if(s=="RMSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RMSE),]
    top.model <- as.character(Summary.Table$Models[which.min(rmse.val)])
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters sorted by Root Mean Square Error (RMSE)")
  }else if(s=="MAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MAE),]
    top.model <- as.character(Summary.Table$Models[which.min(mae.val)])
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters sorted by Mean Absolute Error (MAE)")
  }else if(s=="MSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MSE),]
    top.model <- as.character(Summary.Table$Models[which.min(mse.val)])
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters sorted by Mean Squared Error (MSE)")
  }else if(s=="RAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RAE),]
    top.model <- as.character(Summary.Table$Models[which.min(rae.val)])
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters sorted by Relative Absolute Error (RAE)")
  }else if(s=="AIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$AIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PAIC.val)])
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters sorted by Akaike Information Criterion (AIC)")
  }else if(s=="BIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$BIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PBIC.val)])
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters sorted by Bayesian Information Criterion (BIC)")
  }else if(s=="R2"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$R2,decreasing=TRUE),]
    top.model <- as.character(Summary.Table$Models[which.max(R2.val)])
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters sorted by Coefficient of Determination (R2)")
  }else if(s=="SE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$SE),]
    top.model <- as.character(Summary.Table$Models[which.min(SE.val)])
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters sorted by Standard Error (SE)")
  }else{
    Sort.Summary.Table <- Summary.Table
    message("Summary of Linear External Diffusion Models: List of Error Values and Parameters (Unsorted)")
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
