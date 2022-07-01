#' @title Linear Boyd Intraparticle Diffusion Model
#' @description The Boyd Intraparticle Diffusion Model is frequently applied to analyze if intraparticle diffusion governs the experimental kinetic data. This model assumes that the boundary layer surrounding the adsorbent has a greater effect on the diffusion of solute (Viegas, Campinas, Costa, and Rosa, 2014).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qinf the numerical value for the amount adsorbed at infinite time
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @import utils
#' @return the linear regression and the parameter estimation for the Boyd Intraparticle Diffusion model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <-c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qinf <- 4.8
#' bid.lm(t,qt,qinf)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Boyd, G. E., Adamson, A. W., & Myers, L. S. (1947) <doi:10.1021/ja01203a066> The Exchange Adsorption of Ions from Aqueous Solutions by Organic Zeolites. II. Kinetics1. Journal of the American Chemical Society, 69(11), 2836-2848.
#' @references Viegas, R. M. C., Campinas, M., Costa, H., &; Rosa, M. J. (2014) <doi:10.1007/s10450-014-9617-9> How do the HSDM and Boyd's model compare for estimating intraparticle diffusion coefficients in adsorption processes. Adsorption, 20(5-6), 737-746.
#' @export
bid.lm<- function(t,qt,qinf){
  x <- t
  y <- qt
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(x,y){
    fxnboyd <- y ~ qinf*(1-((6/(pi^2))*((exp(-B*x))+(exp(-B*4*x)/4)+(exp(-B*9*x)/9)+(exp(-B*16*x)/16)+(exp(-B*25*x)/25)
                                        +(exp(-B*36*x)/36)+(exp(-B*49*x)/49)+(exp(-B*64*x)/64)+(exp(-B*81*x)/81)
                                        +(exp(-B*100*x)/100))))
    grdboyd <- data.frame(B = c(0,1),
                          qinf = c(0,1000))
    fit12 <- nls2(fxnboyd,
                  data = dat,
                  start = grdboyd,
                  algorithm = "plinear-random",
                  control = list(maxiter = 1000))
    parsboyd <- as.vector(coefficients(fit12))
    pars_B <- parsboyd[1L]; pars_qinf <- parsboyd[2L];pars_lin <- parsboyd[3L]
    Bmin <- pars_B*0.9; Bmax <- pars_B*1.1
    asymptote <- function(x){ pars_qinf*pars_lin*(1-((6/(pi^2))*((exp(-pars_B*x))+(exp(-pars_B*4*x)/4)+(exp(-pars_B*9*x)/9)+(exp(-pars_B*16*x)/16)+(exp(-pars_B*25*x)/25)
                                                                 +(exp(-pars_B*36*x)/36)+(exp(-pars_B*49*x)/49)+(exp(-pars_B*64*x)/64)+(exp(-pars_B*81*x)/81)
                                                                 +(exp(-pars_B*100*x)/100))))
    }
    qinfmin  <- asymptote(100000) ; qinfmax <- asymptote(100000)*1.1
    Bmin     <- pars_B*0.9        ; Bmax    <- pars_B*1.1
    grdboyd1 <- data.frame(B=c(Bmin,Bmax),
                           qinf=c(qinfmin,qinfmax))
    fit12    <- nls2(fxnboyd,
                     start = grdboyd1,
                     algorithm = "brute-force",
                     control=list(maxiter=1000))
    parsboyd <- as.vector(coefficients(fit12))
    pars_B   <- parsboyd[1L] ;qinf <- parsboyd[2L]
    Bmin     <- pars_B*0.9   ;Bmax <- pars_B*1.1
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
    }
    if(length(Short.F)==0){
      cc<- capture.output(type="message",
                          fit13 <- try(nls2(long.time,
                                                  data=long.data,
                                                  start = grd.long,
                                                  algorithm = "port",
                                                  control=list(maxiter=1000)),
                                       silent=TRUE))
      B.l <- as.vector(coefficients(fit13))
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
      predval <- function(n.fin){
        Col1     <- c(rep(" |",each = n.fin))
        Col2     <- c(rep("|",each = n.fin))
        Bt.val   <- as.vector(B.l*time.l)
        P.Table  <- data.frame(Col1,time.val,Col1,pred.val,Col1,Bt.val,Col2)
        colnames(P.Table) <- c(" |","Time "," |","Pred Val"," |","Bt Val","|")
        message("Estimated Values:")
        print(P.Table, right=T, row.names = F)
      }
    }else{
      cc<- capture.output(type="message",
                          fit14 <- try(nls2(short.time,
                                                  data=short.data,
                                                  start = grd.short,
                                                  algorithm = "port",
                                                  control=list(maxiter=1000)),
                                       silent=TRUE))
      B.s <- as.vector(coefficients(fit14))
      if(is.null(B.s)==TRUE){
        B.s <- ((2*pi)-(((pi^2)*Short.F)/3)-((2*pi)*sqrt(1-(pi*Short.F/3)))/time.s)
      }
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
        predval <- function(n.fin){
          Col1     <- c(rep(" |",each = n.fin))
          Col2     <- c(rep("|",each = n.fin))
          Bt.val   <- as.vector(B.s*time.s)
          P.Table  <- data.frame(Col1,time.val,Col1,pred.val,Col1,Bt.val,Col2)
          colnames(P.Table) <- c(" |","Time "," |","Pred Val"," |","Bt Val","|")
          message("Estimated Values:")
          print(P.Table, right=T, row.names = F)
        }
      }else{
        cc<- capture.output(type="message",
                            fit13 <- try(nls2(long.time,
                                                    data=long.data,
                                                    start = grd.long,
                                                    algorithm = "plinear-random",
                                                    control=list(maxiter=1000)),
                                         silent=TRUE))
        parsboyd.l <- as.vector(coefficients(fit13))
        pars_B <- parsboyd.l[1L]; pars_lin <- parsboyd.l[2L]
        if(is.null(pars_B)==TRUE){
          B.l <- ((-1.4977 + Long.F)/time.l)
        }else{B.l <- pars_B/pars_lin}
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
        predval <- function(n.fin){
          Col1     <- c(rep(" |",each = n.fin))
          Col2     <- c(rep("|",each = n.fin))
          Bt.val   <- c(B.s*time.s,B.l*time.l)
          P.Table  <- data.frame(Col1,time.val,Col1,pred.val,Col1,Bt.val,Col2)
          colnames(P.Table) <- c(" |","Time "," |","Pred Val"," |","Bt Val","|")
          message("Estimated Values")
          print(P.Table, right=T, row.names = F)
        }
      }
    }
    F.val.1<-c(Short.F,Long.F)
    qt.val <-c(F.val.1*qinf)
    if(exists("B.s")==FALSE){
      B.s <-"NA"
    }
    if(exists("B.l")==FALSE){
      B.l <-"NA"
    }
    error <- function(qt.val){
      rmse   <- round((as.numeric(rmse(qt.val,pred.val))),digits=10)
      mae    <- round((as.numeric(mae(qt.val,pred.val))),digits=10)
      mse    <- round((as.numeric(mse(qt.val,pred.val))),digits=10)
      rae    <- round((as.numeric(rae(qt.val,pred.val))),digits=10)
      PAIC   <- round((as.numeric(AIC(BID.lm))),digits=10)
      PBIC   <- round((as.numeric(BIC(BID.lm))),digits=10)
      SE     <- round((as.numeric(sqrt((sum((qt.val-pred.val)^2))/(n.fin-2)))),digits=10)
      rsqtot <- round((as.numeric(summary(BID.lm)$r.squared)),digits=10)
      Col1   <- c(" |"," |"," |"," |"," |"," |"," |"," |")
      Col2   <- c("|","|","|","|","|","|","|","|")
      E.P    <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ","Coefficient of Determination (R2)")
      E.V    <- c(rmse,mae,mse,rae,PAIC,PBIC,SE,rsqtot)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation:")
      print(E.Table, right=F, row.names = F)
    }
    params <- function(qinf){
      param.name <- c("B at short time=","B at long time=")
      param.val <- c(B.s,B.l)
      param.table <- data.frame(param.name,param.val)
      colnames(param.table) <- c("qinf=",qinf)
      message("Boyd Internal Diffusion Parameters")
      print(param.table, right=TRUE, row.names = F)
    }
    message("Boyd Intraparticle Diffusion Model")
    print(summary(BID.lm))
    params(qinf)
    predval(n.fin)
    error(qt.val)
    xval <- seq(min(time.val),max(time.val),length=100)
    eqm  <- slp*xval+ int
    plot.dat <- data.frame (xval,eqm)
    theme_set(theme_bw())
    plot <- ggplot(Bt.table, aes(x=time,y=Bt))+
      geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red")+
      geom_point(size=2)+
      labs(subtitle="Plot of Bt vs time with Boyd intraparticle diffusion model",
           y="Bt",
           x="time",
           title="Boyd Intraparticle Diffusion",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(x,y,qinf){
    fxnboyd <- y ~ qinf*(1-((6/(pi^2))*((exp(-B*x))+(exp(-B*4*x)/4)+(exp(-B*9*x)/9)+(exp(-B*16*x)/16)+(exp(-B*25*x)/25)
                                        +(exp(-B*36*x)/36)+(exp(-B*49*x)/49)+(exp(-B*64*x)/64)+(exp(-B*81*x)/81)
                                        +(exp(-B*100*x)/100))))
    grdboyd <- data.frame(B = c(0,1))
    fit12 <- nls2(fxnboyd,
                  data = dat,
                  start = grdboyd,
                  algorithm = "plinear-random",
                  control = list(maxiter = 1000))
    parsboyd <- as.vector(coefficients(fit12))
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
                          fit13 <- try(nls2(long.time,
                                                  data=long.data,
                                                  start = grd.long,
                                                  algorithm = "port",
                                                  control=list(maxiter=1000)),
                                       silent=TRUE))
      B.l <- as.vector(coefficients(fit13))
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
      predval <- function(n.fin){
        Col1     <- c(rep(" |",each = n.fin))
        Col2     <- c(rep("|",each = n.fin))
        Bt.val   <- as.vector(B.l*time.l)
        P.Table  <- data.frame(Col1,time.val,Col1,pred.val,Col1,Bt.val,Col2)
        colnames(P.Table) <- c(" |","Time "," |","Pred Val"," |","Bt Val","|")
        message("Estimated Values:")
        print(P.Table, right=T, row.names = F)
      }
    }else{
      cc<- utils::capture.output(type="message",
                          fit14 <- try(nls2(short.time,
                                                  data=short.data,
                                                  start = grd.short,
                                                  algorithm = "port",
                                                  control=list(maxiter=1000)),
                                       silent=TRUE))
      B.s <- as.vector(coefficients(fit14))
      if(is.null(B.s)==TRUE){
        B.s <- ((2*pi)-(((pi^2)*Short.F)/3)-((2*pi)*sqrt(1-(pi*Short.F/3)))/time.s)
      }
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
        predval <- function(n.fin){
          Col1     <- c(rep(" |",each = n.fin))
          Col2     <- c(rep("|",each = n.fin))
          Bt.val   <- as.vector(B.s*time.s)
          P.Table  <- data.frame(Col1,time.val,Col1,pred.val,Col1,Bt.val,Col2)
          colnames(P.Table) <- c(" |","Time "," |","Pred Val"," |","Bt Val","|")
          message("Estimated Values:")
          print(P.Table, right=T, row.names = F)
        }
      }else{
        cc<- capture.output(type="message",
                            fit13 <- try(nls2(long.time,
                                                    data=long.data,
                                                    start = grd.long,
                                                    algorithm = "plinear-random",
                                                    control=list(maxiter=1000)),
                                         silent=TRUE))

        parsboyd.l <- as.vector(coefficients(fit13))
        pars_B <- parsboyd.l[1L]; pars_lin <- parsboyd.l[2L]
        if(is.null(pars_B)==TRUE){
          B.l <- ((-1.4977 + Long.F)/time.l)
        }else{B.l <- pars_B/pars_lin}
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
        predval <- function(n.fin){
          Col1     <- c(rep(" |",each = n.fin))
          Col2     <- c(rep("|",each = n.fin))
          Bt.val   <- c(B.s*time.s,B.l*time.l)
          P.Table  <- data.frame(Col1,time.val,Col1,pred.val,Col1,Bt.val,Col2)
          colnames(P.Table) <- c(" |","Time "," |","Pred Val"," |","Bt Val","|")
          message("Estimated Values:")
          print(P.Table, right=T, row.names = F)
        }
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
    error <- function(qt.val){
      rmse   <- round((as.numeric(rmse(qt.val,pred.val))),digits=10)
      mae    <- round((as.numeric(mae(qt.val,pred.val))),digits=10)
      mse    <- round((as.numeric(mse(qt.val,pred.val))),digits=10)
      rae    <- round((as.numeric(rae(qt.val,pred.val))),digits=10)
      PAIC   <- round((as.numeric(AIC(BID.lm))),digits=10)
      PBIC   <- round((as.numeric(BIC(BID.lm))),digits=10)
      SE     <- round((as.numeric(sqrt((sum((qt.val-pred.val)^2))/(n.fin-2)))),digits=10)
      rsqtot <- round((as.numeric(summary(BID.lm)$r.squared)),digits=10)
      Col1   <- c(" |"," |"," |"," |"," |"," |"," |"," |")
      Col2   <- c("|","|","|","|","|","|","|","|")
      E.P    <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ","Coefficient of Determination (R2)")
      E.V    <- c(rmse,mae,mse,rae,PAIC,PBIC,SE,rsqtot)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation:")
      print(E.Table, right=F, row.names = F)
    }
    params <- function(qinf){
      param.name <- c("B at short time=","B at long time=")
      param.val <- c(B.s,B.l)
      param.table <- data.frame(param.name,param.val)
      colnames(param.table) <- c("qinf=",qinf)
      message("Boyd Internal Diffusion Parameters")
      print(param.table, right=TRUE, row.names = F)
    }
    message("Boyd Intraparticle Diffusion Model")
    print(summary(BID.lm))
    params(qinf)
    predval(n.fin)
    error(qt.val)
    xval <- seq(min(time.val),max(time.val),length=100)
    eqm  <- slp*xval+ int
    plot.dat <- data.frame (xval,eqm)
    theme_set(theme_bw())
    plot <- ggplot(Bt.table, aes(x=time,y=Bt))+
      geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red")+
      geom_point(size=2)+
      labs(subtitle="Plot of Bt vs time with Boyd intraparticle diffusion model",
           y="Bt",
           x="time",
           title="Boyd Intraparticle Diffusion",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  if(missing(qinf)){
    EQ1(x,y)}
  else if(is.null(qinf)){
    EQ1(x,y)}
  else if(isFALSE(qinf)){
    EQ1(x,y)}
  else{EQ2(x,y,qinf)}
}
