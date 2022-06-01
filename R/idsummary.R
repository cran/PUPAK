#' @title Internal Diffusion Models Summary
#' @description Summarized results of parameter and error values collected from internal diffusion models, namely: Boyd Internal Diffusion, Crank, and Weber and Morris
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qinf the numerical value for the amount adsorbed at infinite time
#' @param sort.by the name of the statistical error parameter in which the models are sorted in either increasing or decreasing order. The only accepted arguments are "RMSE" for Relative Mean Square Error, 'MAE' for Mean Absolute Error, 'MSE' for Mean Squared Error, 'RAE' for Relative Absolute Error, 'AIC' for Akaike Information Criterion, 'BIC' for Bayesian Information Criterion, 'R2' for Coefficient of Determination, and 'SE' for Standard Error Estimate. This argument is case-sensitive, and failure to input the correct value will yield a summary of models in alphabetical order.
#' @import nls2
#' @import stats
#' @import Metrics
#' @import segmented
#' @import utils
#' @return the summarized error and parameter values from internal diffusion models.
#' @examples
#' \donttest{
#' t  <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' idsummary(t,qt,qinf=4.8,"SE")}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @export
idsummary <- function(t,qt,qinf,sort.by){
  errors <- new.env()
  parameters <- new.env()
  summary.data <- new.env()
  linearfunction <- new.env()
  if(missing(qinf)){
    stop("qinf is required to run the function. The qe value can be solve using isotherm models.
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
  t <- t; qt <- qt;  qinf <- qinf
  Crank.sum.nl <- function(t,qt,qinf){
    x <- t ;y <- qt ;qinf <- qinf
    dat  <- data.frame(x,y)
    n.dat  <- nrow(na.omit(dat))
    fxncrnk <- y ~ qinf*(1-((6/(pi^2))*((exp(-Dc*(pi^2)*x/(r^2)))+((exp(-Dc*4*(pi^2)*x/(r^2)))/4)+((exp(-Dc*9*(pi^2)*x/(r^2)))/9)
                                        +((exp(-Dc*16*(pi^2)*x/(r^2)))/16)+((exp(-Dc*25*(pi^2)*x/(r^2)))/25)+((exp(-Dc*36*(pi^2)*x/(r^2)))/36)
                                        +((exp(-Dc*49*(pi^2)*x/(r^2)))/49)+((exp(-Dc*64*(pi^2)*x/(r^2)))/64)+((exp(-Dc*81*(pi^2)*x/(r^2)))/81)
                                        +((exp(-Dc*100*(pi^2)*x/(r^2)))/100))))
    grdcrnk <- data.frame(Dc = c(0,10),
                          r=c(0,10))
    fit33 <- nls2(fxncrnk,
                  data = dat,
                  start = grdcrnk,
                  algorithm = "plinear-random",
                  control = list(maxiter = 2000))
    parscrnk <- as.vector(coefficients(fit33))
    pars_Dc <- parscrnk[1L]; pars_r <- parscrnk[2L]; pars_lin <- parscrnk[3L]
    Dcmin <- pars_Dc*0.9 ; Dcmax <- pars_Dc*1.1
    rmin <- pars_r*0.9   ; rmax <- pars_r*1.1
    grdcrnk1 <- data.frame(Dc=c(Dcmin,Dcmax),
                           r=c(rmin,rmax))
    fit33 <- nls2(fxncrnk,
                  start = grdcrnk1,
                  algorithm = "brute-force",
                  control=list(maxiter=1000))
    summary.data$Crank<-summary(fit33)
    errors$rmse.Crank <- rmse(y,predict(fit33))
    errors$mae.Crank  <- mae(y,predict(fit33))
    errors$mse.Crank  <- mse(y,predict(fit33))
    errors$rae.Crank  <- rae(y,predict(fit33))
    errors$PAIC.Crank <- AIC(fit33)
    errors$PBIC.Crank <- BIC(fit33)
    errors$SE.Crank   <- sqrt((sum((y-predict(fit33))^2))/(n.dat-2))
    parscrnk <- as.vector(coefficients(fit33))
    parameters$Crank.Dc <- parscrnk[1L]
    parameters$Crank.r <- parscrnk[2L]
  }
  WAM.sum.pl<- function(t,qt){
    x <- sqrt(t) ;y <- qt
    dat<- data.frame(x,y)
    n  <- nrow(na.omit(dat))
    WM.lm <- lm(y ~ x)
    fit34 <- as.character(segmented(WM.lm, seg.Z=~x, psi=x[3], control=seg.control(display=FALSE,)), K=3)
    if(length(fit34)==12){
      fit35 <- lm(y ~ x)
      lin.val <- coef(fit35)
      int<-as.numeric(lin.val[1])
      slp <-as.numeric(lin.val[2])
      pred.val <- ((slp*x)+int)
      summary.data$Weber.and.Morris<-summary(fit35)
      parameters$WAM.kWAM1 <- slp
      parameters$WAM.C1 <- int
      errors$rmse.WAM <- rmse(y,pred.val)
      errors$mae.WAM  <- mae(y,pred.val)
      errors$mse.WAM  <- mse(y,pred.val)
      errors$rae.WAM  <- rae(y,pred.val)
      errors$PAIC.WAM <- AIC(fit35)
      errors$PBIC.WAM <- BIC(fit35)
      errors$SE.WAM   <- sqrt((sum((y-pred.val)^2))/(n-2))
      errors$R2.WAM   <- cor(x,y)^2
      brk.pts <- NULL
    }else{
      fit34 <- segmented(WM.lm, seg.Z=~x, psi=x[3], control=seg.control(display=FALSE,), K=3)
      brk.pts <- fit34$psi[, 2]
      seg.range <- fit34$rangeZ[, 1]
      b <- intercept(fit34); y.int <- unlist(b)
      a <- slope(fit34)    ; slp <- unlist(a)
      int1 <- as.numeric(y.int[1])
      slp1 <- as.numeric(slp[1])
      int2 <- as.numeric(y.int[2])
      slp2 <- as.numeric(slp[2])
      brk <- brk.pts[1]
      firstlinearx <- x[which(x < brk )]
      secondlinearx <- x[which(x > brk)]
      firstlineary <- y[which(dat$x < brk )]
      secondlineary <- y[which(dat$x > brk)]
      low.pred <- ((slp1*firstlinearx)+int1)
      high.pred <- ((slp2*secondlinearx)+int2)
      pred.val <- c(low.pred,high.pred)
      RSS1 <- sum((firstlineary-low.pred)^2)
      RSS2 <- sum((secondlineary-high.pred)^2)
      TSS1 <- sum((firstlineary-mean(firstlineary))^2)
      TSS2 <- sum((secondlineary-mean(secondlineary))^2)
      Breakpoint <- function(brk.pts){
        capture.output(message("x=",brk.pts),type="message")
      }
      summary.data$Weber.and.Morris<-summary(fit34)
      parameters$WAM.kWAM1 <- slp1
      parameters$WAM.kWAM2 <- slp2
      parameters$WAM.C1 <- int1
      parameters$WAM.C2 <- int2
      errors$rmse.WAM <- rmse(y,pred.val)
      errors$mae.WAM  <- mae(y,pred.val)
      errors$mse.WAM  <- mse(y,pred.val)
      errors$rae.WAM  <- rae(y,pred.val)
      errors$PAIC.WAM <- AIC(fit34)
      errors$PBIC.WAM <- BIC(fit34)
      errors$SE.WAM   <- sqrt((sum((y-pred.val)^2))/(n-2))
      errors$R21.WAM  <- as.numeric(cor(firstlinearx,firstlineary)^2)
      errors$R22.WAM  <- as.numeric(cor(secondlinearx,secondlineary)^2)
      errors$R2.WAM   <- ((1-((RSS1 + RSS2)/(TSS1 + TSS2))))
      linearfunction$Weber.and.Morris.breakpoint <- brk.pts
    }
  }
  BID.sum.l<- function(t,qt,qinf){
    x <- t ;y <- qt
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    fxnboyd <- y ~ qinf*(1-((6/(pi^2))*((exp(-B*x))+(exp(-B*4*x)/4)+(exp(-B*9*x)/9)+(exp(-B*16*x)/16)+(exp(-B*25*x)/25)
                                        +(exp(-B*36*x)/36)+(exp(-B*49*x)/49)+(exp(-B*64*x)/64)+(exp(-B*81*x)/81)
                                        +(exp(-B*100*x)/100))))
    grdboyd <- data.frame(B = c(0,1))
    fit36 <- nls2(fxnboyd,
                  data = dat,
                  start = grdboyd,
                  algorithm = "plinear-random",
                  control = list(maxiter = 1000))
    parsboyd <- as.vector(coefficients(fit36))
    pars_B <- parsboyd[1L]
    Bmin <- pars_B*0.9; Bmax <- pars_B*1.1
    F.val   <- as.vector(y/qinf)
    Short.F <- as.vector(F.val[which(F.val <= 0.85 )])
    Long.F  <- as.vector(F.val[which(F.val >= 0.85 )])
    Long.F  <- as.vector(Long.F[which(Long.F <= 1)])
    time.s  <- as.vector(x[which(F.val <= 0.85)])
    time.l  <- as.vector(x[which(F.val > 0.85)])
    time.l  <- as.vector(time.l[which(Long.F <= 1)])
    short.data <- data.frame(Short.F,time.s)
    short.time <- time.s~(((2*pi)-(((pi^2)*Short.F)/3)-((2*pi)*(1-(pi*Short.F/3))^(0.5)))/B.short)
    grd.short  <- data.frame(B.short=c(Bmin,Bmax))
    long.data  <- data.frame(Long.F,time.l)
    long.time  <- time.l~((-log(((pi^2)/6)*(1-Long.F)))/B.long)
    grd.long   <- data.frame(B.long=c(Bmin,Bmax))
    if(sum(short.data$Short.F)==0){
      F.val   <- F.val[which(F.val > 0)]
      Short.F <- as.vector(F.val[which(F.val <= 0.85 )])
      time.s  <- x[which(F.val <= 0.85)]
      short.data <- data.frame(Short.F,time.s)
    }else{}
    if(length(Short.F)==0){
      cc<- capture.output(type="message",
                          fit37 <- try(nls2::nls2(long.time,
                                                  data=long.data,
                                                  start = grd.long,
                                                  algorithm = "port",
                                                  control=list(maxiter=1000)),
                                       silent=TRUE))
      B.l <- as.vector(coefficients(fit37))
      exp.Bt.long <- as.vector(-log(((pi^2)/6)*(1-Long.F)))
      Bt.table <- data.frame(exp.Bt.long,time.l)
      colnames(Bt.table) <- c("Bt","time")
      time <- Bt.table$time
      Bt   <- Bt.table$Bt
      n.fin  <- nrow((Bt.table))
      BID.lm <- lm(Bt~time)
      lm.val <- as.vector(coefficients(BID.lm))
      slp <- lm.val[2L]; int <-lm.val[1L]
      qt.pred.l  <- function(B,t,qinf){qinf*(1-(((pi^2)/6)*exp(-B*t)))}
      pred.val <- qt.pred.l(B.l,time.l,qinf)
      time.val <- time.l
    }else{
      cc<- capture.output(type="message",
                          fit38 <- try(nls2::nls2(short.time,
                                                  data=short.data,
                                                  start = grd.short,
                                                  algorithm = "port",
                                                  control=list(maxiter=1000)),
                                       silent=TRUE))
      B.s <- as.vector(coefficients(fit38))
      if(length(Long.F)==0){
        exp.Bt.short <- as.vector(((2*pi)-(((pi^2)*Short.F)/3)-((2*pi)*(1-(pi*Short.F/3))^(0.5))))
        Bt.table <- data.frame(exp.Bt.short,time.s)
        colnames(Bt.table) <- c("Bt","time")
        time <- Bt.table$time
        Bt   <- Bt.table$Bt
        n.fin  <- nrow(na.omit(Bt.table))
        BID.lm <- lm(Bt~time)
        lm.val <- as.vector(coefficients(BID.lm))
        slp <- lm.val[2L]; int <-lm.val[1L]
        qt.pred.s <- function(B,t,qinf){(((6/pi^(3/2))*sqrt(B*t))-((3/(pi^2))*B*t))*qinf}
        pred.val <- qt.pred.s(B.s,time.s,qinf)
        time.val <- time.s
      }else{
        cc<- capture.output(type="message",
                            fit37 <- try(nls2::nls2(long.time,
                                                    data=long.data,
                                                    start = grd.long,
                                                    algorithm = "plinear-random",
                                                    control=list(maxiter=1000)),
                                         silent=TRUE))

        parsboyd.l <- as.vector(coefficients(fit37))
        pars_B <- parsboyd.l[1L]; pars_lin <- parsboyd.l[2L]
        B.l <- pars_B/pars_lin
        time.l <- long.data$time.l
        exp.Bt.short <- as.vector(((2*pi)-(((pi^2)*Short.F)/3)-((2*pi)*(1-(pi*Short.F/3))^(0.5))))
        exp.Bt.long  <- as.vector((-log(((pi^2)/6)*(1-Long.F))))
        Bt.s   <- data.frame(exp.Bt.short,time.s)
        Bt.l   <- data.frame(exp.Bt.long,time.l)
        colnames(Bt.s) <- c("Bt","time")
        colnames(Bt.l) <- c("Bt","time")
        Bt.table <- rbind(Bt.s,Bt.l)
        time   <- Bt.table$time
        Bt     <- Bt.table$Bt
        n.fin  <- nrow(na.omit(Bt.table))
        BID.lm <- lm(Bt~time)
        lm.val <- as.vector(coefficients(BID.lm))
        slp <- lm.val[2L]; int <-lm.val[1L]
        qt.pred.s <- function(B,t,qinf){(((6/pi^(3/2))*sqrt(B*t))-((3/(pi^2))*B*t))*qinf}
        qt.pred.l <- function(B,t,qinf){qinf*(1-(((pi^2)/6)*exp(-B*t)))}
        qt.s     <- qt.pred.s(B.s,time.s,qinf)
        qt.l     <- qt.pred.l(B.l,time.l,qinf)
        pred.val <- c(qt.s,qt.l)
        time.val <- c(time.s,time.l)
      }
    }
    F.val.1<-c(Short.F,Long.F)
    qt.val <-c(F.val.1*qinf)
    if(exists("B.s")==FALSE){
      B.s <-"NA"
    }else{}
    if(exists("B.l")==FALSE){
      B.l <-"NA"
    }else{}
    summary.data$Boyd.Intraparticle.Diffusion<-summary(BID.lm)
    parameters$BID.B.short <- B.s
    parameters$BID.B.long <- B.l
    errors$rmse.BID <- rmse(qt.val,pred.val)
    errors$mae.BID  <- mae(qt.val,pred.val)
    errors$mse.BID  <- mse(qt.val,pred.val)
    errors$rae.BID  <- rae(qt.val,pred.val)
    errors$PAIC.BID <- AIC(BID.lm)
    errors$PBIC.BID <- BIC(BID.lm)
    errors$SE.BID   <- sqrt((sum((qt.val-pred.val)^2))/(n.fin-2))
    errors$R2.BID   <- cor(time,Bt)^2
  }

  Crank.sum.nl(t,qt,qinf)
  WAM.sum.pl(t,qt)
  BID.sum.l(t,qt,qinf)
  Crank.Dc    <- format(as.numeric(parameters[["Crank.Dc"]]),scientific = 5)
  Crank.r     <- format(as.numeric(parameters[["Crank.r"]]),scientific = 5)
  WAM.kWAM1   <- format(as.numeric(parameters[["WAM.kWAM1"]]),scientific = 5)
  WAM.C1      <- format(as.numeric(parameters[["WAM.C1"]]),scientific = 5)
  if(parameters[["BID.B.long"]]=="NA"){
    BID.B.long  <- "-"
  }else{
    BID.B.long  <- format(as.numeric(parameters[["BID.B.long"]]),scientific = 5)
  }
  if(parameters[["BID.B.short"]]=="NA"){
    BID.B.short <- "-"
  }else{
    BID.B.short <- format(as.numeric(parameters[["BID.B.short"]]),scientific = 5)
  }
  if(length(parameters[["WAM.kWAM2"]])==0){
    WAM.kWAM2   <- "-"
  }else{
    WAM.kWAM2   <- format(as.numeric(parameters[["WAM.kWAM2"]]),scientific = 5)
  }
  if(length(parameters[["WAM.C2"]])==0){
    WAM.C2 <- "-"
  }else{
    WAM.C2      <- format(as.numeric(parameters[["WAM.C2"]]),scientific = 5)
  }
  param1.name <- c("B.short.time","Dc","kWAM1")
  param2.name <- c("B.long.time","r","C1")
  param3.name <- c("-","-","KwaM2")
  param4.name <- c("-","-","C2")
  param1.val <- as.character(c(BID.B.short,Crank.Dc,WAM.kWAM1))
  param2.val <- as.character(c(BID.B.long,Crank.r,WAM.C1))
  param3.val <- as.character(c("-","-",WAM.kWAM2))
  param4.val <- as.character(c("-","-",WAM.C2))
  R2.BID <- round(as.numeric(errors[["R2.BID"]]),digits = 5)
  R2.WAM <- round(as.numeric(errors[["R2.WAM"]]),digits = 5)
  Models <- c("BID","Crank","WAM")
  rmse.val <- as.numeric(c(errors[["rmse.BID"]],errors[["rmse.Crank"]],errors[["rmse.WAM"]]))
  mae.val  <- as.numeric(c(errors[["mae.BID"]],errors[["mae.Crank"]],errors[["mae.WAM"]]))
  mse.val  <- as.numeric(c(errors[["mse.BID"]],errors[["mse.Crank"]],errors[["mse.WAM"]]))
  rae.val  <- as.numeric(c(errors[["rae.BID"]],errors[["rae.Crank"]],errors[["rae.WAM"]]))
  PAIC.val <- as.numeric(c(errors[["PAIC.BID"]],errors[["PAIC.Crank"]],errors[["PAIC.WAM"]]))
  PBIC.val <- as.numeric(c(errors[["PBIC.BID"]],errors[["PBIC.Crank"]],errors[["PBIC.WAM"]]))
  SE.val   <- as.numeric(c(errors[["SE.BID"]],errors[["SE.Crank"]],errors[["SE.WAM"]]))
  R2.val   <- as.character(c(R2.BID,"-",R2.WAM))
  rmse.val <- as.character(round(rmse.val,digits=5))
  mae.val  <- as.character(round(mae.val, digits=5))
  mse.val  <- as.character(round(mse.val, digits=5))
  rae.val  <- as.character(round(rae.val, digits=5))
  PAIC.val <- as.character(round(PAIC.val,digits=4))
  PBIC.val <- as.character(round(PBIC.val,digits=4))
  SE.val   <- as.character(round(SE.val,  digits=5))
  Summary.Table <- data.frame(Models,rmse.val,mae.val,mse.val,rae.val,PAIC.val,PBIC.val,R2.val,SE.val,param1.name,param1.val,param2.name,param2.val,param3.name,param3.val,param4.name,param4.val)
  colnames(Summary.Table) <- c("Models","RMSE","MAE","MSE","RAE","AIC","BIC","R2","SE","Parameter1","       ","Parameter2","        ","Parameter 3","          ","Parameter 4","      ")
  if(s=="RMSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RMSE),]
    top.model <- as.character(Summary.Table$Models[which.min(rmse.val)])
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters sorted by Root Mean Square Error (RMSE)")
  }else if(s=="MAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MAE),]
    top.model <- as.character(Summary.Table$Models[which.min(mae.val)])
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters sorted by Mean Absolute Error (MAE)")
  }else if(s=="MSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MSE),]
    top.model <- as.character(Summary.Table$Models[which.min(mse.val)])
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters sorted by Mean Squared Error (MSE)")
  }else if(s=="RAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RAE),]
    top.model <- as.character(Summary.Table$Models[which.min(rae.val)])
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters sorted by Relative Absolute Error (RAE)")
  }else if(s=="AIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$AIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PAIC.val)])
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters sorted by Akaike Information Criterion (AIC)")
  }else if(s=="BIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$BIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PBIC.val)])
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters sorted by Bayesian Information Criterion (BIC)")
  }else if(s=="R2"){
    Summary.Table[2,8] <- -1
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$R2,decreasing=TRUE),]
    top.model <- as.character(Summary.Table$Models[which.max(Summary.Table$R2)])
    Sort.Summary.Table[3,8] <- "-"
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters sorted by Coefficient of Determination (R2)")
  }else if(s=="SE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$SE),]
    top.model <- as.character(Summary.Table$Models[which.min(SE.val)])
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters sorted by Standard Error (SE)")
  }else{
    Sort.Summary.Table <- Summary.Table
    message("Summary of Internal Diffusion Models: List of Error Values and Parameters (Unsorted)")
    top.model <- as.character("Unsorted")
  }
  Param.Table <- Sort.Summary.Table[c(-2,-3,-4,-5,-6,-7,-8,-9)]
  Sort.Summary.Table <-Sort.Summary.Table[c(-10,-11,-12,-13,-14,-15,-16,-17)]
  print(Sort.Summary.Table, right=FALSE,row.names = FALSE)
  print(Param.Table, right=FALSE,row.names = FALSE )
  if(as.character(top.model)!= "Unsorted"){
    if(top.model=="BID"){
      top.model <- "Boyd.Intraparticle.Diffusion"
    }
    if(top.model=="WAM"){
      top.model <- "Weber.and.Morris"
    }
    message("Summary of the most fitted model (",top.model,")",sep="")
    if(top.model!="Crank"){
      if(top.model=="Weber.and.Morris"){
        break.point <-data.frame(c(""),c(""))
        colnames(break.point) <- c("Break Point =",linearfunction[["Weber.and.Morris.breakpoint"]])
        print(break.point, right=TRUE, row.names = F)
        if(length((errors[["R21.WAM"]]))==1){
          CoefR2 <-data.frame(c("Second Equation R2="),c(errors[["R22.WAM"]]))
          colnames(CoefR2) <- c("First Equation  R2=",errors[["R21.WAM"]])
          message("Coefficient of Determination (R2) for Equation 1 and Equation 2")
          print(CoefR2, right=TRUE, row.names = F)
        }
      }
    }
    print(summary.data[[top.model]])
  }else{}
}
