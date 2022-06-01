#' @title Non-Linear External Diffusion Models Summary
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
#' @import nls2
#' @import stats
#' @import Metrics
#' @import utils
#' @return the summarized error and parameter values from external diffusion models.
#' @examples
#' \donttest{
#' t  <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' Ct <- c(10.000,8.141,8.056,7.949,7.863,7.799,7.778,7.756,7.692)
#' edsummary.nl(t,qt,Ct,qinf=4.8,Co=10,m=0.05,V=0.1,b=1.3,"SE")}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @export
edsummary.nl<- function(t,qt,Ct,qinf,Co,m,V,b,sort.by){
  errors <- new.env()
  parameters <- new.env()
  summary.data <- new.env()
  t<-t     ; b <- b
  qt<-qt   ; qinf <-qinf
  Ct <- Ct ; Co <- Co
  m <- m   ; V <- V
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
  if(missing(sort.by)){
    s<-"sort.by"
  }else if(is.null(sort.by)){
    s<-"sort.by"
  }else if(isFALSE(sort.by)){
    s<-"sort.by"
  }else{
    s<-sort.by
  }
  BED.sum.nl <- function(t,qt,qinf){
    x <- t ;y <- qt
    qinf <- qinf
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    fxnbed <- y ~ (qinf*(1 - exp(-R*x)))
    grdbed <- data.frame(R = c(0,10))
    cc<- capture.output(type="message",
                        fit30 <- try(nls2::nls2(fxnbed,
                                                data = dat,
                                                start = grdbed,
                                                algorithm = "port",
                                                control = list(maxiter = 1000)),
                                     silent=TRUE))
    if(is.null(fit30)==TRUE){
      fit30 <- nls2(fxnbed,
                    data = dat,
                    start = grdbed,
                    algorithm = "plinear-random",
                    control = list(maxiter = 1000))
      parsbed <- as.vector(coefficients(fit30))
      pars_R <- parsbed[1L]; pars_lin <- parsbed[2L]
      Rmin <- pars_R*0.9; Rmax <- pars_R*1.1
      grdbed1 <- data.frame(R=c(Rmin,Rmax))
      fit30 <- nls2(fxnbed,
                    start = grdbed1,
                    algorithm = "brute-force",
                    control=list(maxiter=1000))
    }else{}
    summary.data$Boyd.External.Diffusion<-summary(fit30)
    errors$rmse.BED <- rmse(y,predict(fit30))
    errors$mae.BED  <- mae(y,predict(fit30))
    errors$mse.BED  <- mse(y,predict(fit30))
    errors$rae.BED  <- rae(y,predict(fit30))
    errors$PAIC.BED <- AIC(fit30)
    errors$PBIC.BED <- BIC(fit30)
    errors$SE.BED   <- sqrt((sum((y-predict(fit30))^2))/(n.dat-2))
    parsbed <- as.vector(coefficients(fit30))
    parameters$BED.R <- parsbed[1L]
  }
  FAS.sum.nl <- function(t,Ct,m,V,b,Co){
    x <- t        ; y <- Ct
    m <- m; V <- V; b <- b; Co <- Co
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    fxnfas <- y ~Co*((1/(1 + (m/V)*b)) + (((m/V)*b)/(1 + (m/V)*b)) * exp(-((1 + (m/V)*b)/((m/V)*b))*kfasS*x))
    grdfas <- data.frame(kfasS  = c(0,10))
    cc<- capture.output(type="message",
                        fit31 <- try(nls2::nls2(fxnfas,
                                                data = dat,
                                                start = grdfas,
                                                algorithm = "port",
                                                control = list(maxiter = 1000)),
                                     silent=TRUE))
    if(is.null(fit31)==TRUE){
      fit31 <- nls2(fxnfas,
                    data = dat,
                    start = grdfas,
                    algorithm = "plinear-random",
                    control = list(maxiter = 1000))
      parsfas <- as.vector(coefficients(fit31))
      pars_kfasS <- parsfas[1L]
      kfasSmin <-  pars_kfasS*0.9 ;kfasSmax <-  pars_kfasS*1.1
      grdem1 <- data.frame(kfasS = c(kfasSmin,kfasSmax))
      fit31 <- nls2(fxnfas,
                    start = grdem1,
                    algorithm = "brute-force",
                    control=list(maxiter=1000))
    }else{}
    summary.data$Furusawa.and.Smith<-summary(fit31)
    errors$rmse.FAS <- rmse(y,predict(fit31))
    errors$mae.FAS  <- mae(y,predict(fit31))
    errors$mse.FAS  <- mse(y,predict(fit31))
    errors$rae.FAS  <- rae(y,predict(fit31))
    errors$PAIC.FAS <- AIC(fit31)
    errors$PBIC.FAS <- BIC(fit31)
    errors$SE.FAS   <- sqrt((sum((y-predict(fit31))^2))/(n.dat-2))
    parsfas <- as.vector(coefficients(fit31))
    parameters$FAS.kFASS <- parsfas[1L]
  }
  MAW.sum.nl <- function(t,Ct,Co){
    x <- t   ;y <- Ct
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    fxnmaw <- y ~ Co*exp(-kmaw*x)
    grdmaw <- data.frame(kmaw  = c(0,1))
    cc<- capture.output(type="message",
                        fit32 <- try(nls2::nls2(fxnmaw,
                                                data = dat,
                                                start = grdmaw,
                                                algorithm = "port",
                                                control = list(maxiter = 1000)),
                                     silent=TRUE))
    if(is.null(fit32)==TRUE){
      fit32 <- nls2(fxnmaw,
                    data = dat,
                    start = grdmaw,
                    algorithm = "plinear-random",
                    control = list(maxiter = 1000))
      parsmaw <- as.vector(coefficients(fit32))
      pars_kmaw <- parsmaw[1L]
      kmawmin <-  pars_kmaw*0.9 ;kmawmax <-  pars_kmaw*1.1
      grdem1 <- data.frame(kmaw=c(kmawmin,kmawmax))
      fit32 <- nls2(fxnmaw,
                    start = grdem1,
                    algorithm = "brute-force",
                    control=list(maxiter=1000))
    }else{}
    summary.data$Matthews.and.Weber<-summary(fit32)
    errors$rmse.MAW <- rmse(y,predict(fit32))
    errors$mae.MAW  <- mae(y,predict(fit32))
    errors$mse.MAW  <- mse(y,predict(fit32))
    errors$rae.MAW  <- rae(y,predict(fit32))
    errors$PAIC.MAW <- AIC(fit32)
    errors$PBIC.MAW <- BIC(fit32)
    errors$SE.MAW   <- sqrt((sum((y-predict(fit32))^2))/(n.dat-2))
    parsmaw <- as.vector(coefficients(fit32))
    parameters$MAW.kMAWS <- parsmaw[1L]
  }
  BED.sum.nl(t,qt,qinf)
  FAS.sum.nl(t,Ct,m,V,b,Co)
  MAW.sum.nl(t,Ct,Co)
  BED.R <- parameters[["BED.R"]]
  FAS.kFASS <- parameters[["FAS.kFASS"]]
  MAW.kMAWS <- parameters[["BED.R"]]
  param1.name <- c("R","kFASS","kMAWS")
  param1.val <- as.numeric(c(BED.R,FAS.kFASS,MAW.kMAWS))
  param1.val <- round(param1.val,digits=6)
  Models <- c("Boyd.External.Diffusion","Furusawa.and.Smith","Matthews.and.Weber")
  rmse.val <- as.numeric(c(errors[["rmse.BED"]],errors[["rmse.FAS"]],errors[["rmse.MAW"]]))
  mae.val  <- as.numeric(c(errors[["mae.BED"]],errors[["mae.FAS"]],errors[["mae.MAW"]]))
  mse.val  <- as.numeric(c(errors[["mse.BED"]],errors[["mse.FAS"]],errors[["mse.MAW"]]))
  rae.val  <- as.numeric(c(errors[["rae.BED"]],errors[["rae.FAS"]],errors[["rae.MAW"]]))
  PAIC.val <- as.numeric(c(errors[["PAIC.BED"]],errors[["PBIC.FAS"]],errors[["PBIC.MAW"]]))
  PBIC.val <- as.numeric(c(errors[["PAIC.BED"]],errors[["PAIC.FAS"]],errors[["PAIC.MAW"]]))
  SE.val   <- as.numeric(c(errors[["SE.BED"]],errors[["SE.FAS"]],errors[["SE.MAW"]]))
  rmse.val <- round(rmse.val,digits=6)
  mae.val  <- round(mae.val,digits=6)
  mse.val  <- round(mse.val,digits=6)
  rae.val  <- round(rae.val,digits=6)
  PAIC.val <- round(PAIC.val,digits=6)
  PBIC.val <- round(PBIC.val,digits=6)
  SE.val   <- round(SE.val,digits=6)
  Summary.Table <- data.frame(Models,rmse.val,mae.val,mse.val,rae.val,PAIC.val,PBIC.val,SE.val,param1.name,param1.val)
  colnames(Summary.Table) <- c("Models","RMSE","MAE","MSE","RAE","AIC","BIC","SE","Parameter","     ")
  if(s=="RMSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RMSE),]
    top.model <- as.character(Summary.Table$Models[which.min(rmse.val)])
    message("Summary of Non-Linear External Diffusion Models: List of Error Values and Parameters sorted by Root Mean Square Error (RMSE)")
  }else if(s=="MAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MAE),]
    top.model <- as.character(Summary.Table$Models[which.min(mae.val)])
    message("Summary of Non-Linear External Diffusion Models: List of Error Values and Parameters sorted by Mean Absolute Error (MAE)")
  }else if(s=="MSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MSE),]
    top.model <- as.character(Summary.Table$Models[which.min(mse.val)])
    message("Summary of Non-Linear External Diffusion Models: List of Error Values and Parameters sorted by Mean Squared Error (MSE)")
  }else if(s=="RAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RAE),]
    top.model <- as.character(Summary.Table$Models[which.min(rae.val)])
    message("Summary of Non-Linear External Diffusion Models: List of Error Values and Parameters sorted by Relative Absolute Error (RAE)")
  }else if(s=="AIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$AIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PAIC.val)])
    message("Summary of Non-Linear External Diffusion Models: List of Error Values and Parameters sorted by Akaike Information Criterion (AIC)")
  }else if(s=="BIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$BIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PBIC.val)])
    message("Summary of Non-Linear External Diffusion Models: List of Error Values and Parameters sorted by Bayesian Information Criterion (BIC)")
  }else if(s=="SE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$SE),]
    top.model <- as.character(Summary.Table$Models[which.min(SE.val)])
    message("Summary of Non-Linear External Diffusion Models: List of Error Values and Parameters sorted by Standard Error (SE)")
  }else{
    Sort.Summary.Table <- Summary.Table
    message("Summary of Non-Linear External Diffusion Models: List of Error Values and Parameters (Unsorted)")
    top.model <- as.character("Unsorted")
  }
  Param.Table <- Sort.Summary.Table[c(-2,-3,-4,-5,-6,-7,-8)]
  Sort.Summary.Table <-Sort.Summary.Table[c(-9,-10)]
  print(Sort.Summary.Table, right=FALSE,row.names = FALSE)
  print(Param.Table, right=FALSE,row.names = FALSE )
  if(as.character(top.model)!= "Unsorted"){
    message("Summary of the most fitted model(",top.model,")",sep="")
    print(summary.data[[top.model]])
  }
}
