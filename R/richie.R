#' @title Non-Linear Ritchie's Equation
#' @description The Ritchie's Equation is a well-known empirical rate equation for gas particle adsorption on solid surfaces. Assuming that the rate of adsorption is exclusively determined by the percentage of vacant sites at time t (Kaki, Gögsu, Altindal, Salih, and Bekaroglu, 2020).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qe the numerical value for the amount adsorbed at equilibrium. If this parameter is not defined, it will be estimated.
#' @param n the Richie's equation order of reaction. If the parameter value n=1, the function will proceed to first-order Richie's equation. If the value n=2, the function will proceed to second-order Richie's equation, and if the n is not defined, the value of n will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameters estimation for the Richie adsorption kinetic model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qe <- 4.68
#' richie(t,qt,qe)
#' richie(t,qt,n=3)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Ritchie, A. G. (1977) <doi:10.1039/F19777301650> Alternative to the Elovich equation for the kinetics of adsorption of gases on solids. Journal of the Chemical Society, Faraday Transactions 1: Physical Chemistry in Condensed Phases, 73, 1650-1653.
#' @references Kaki, E., Gögsu, N., Altindal, A., Salih, B., &; Bekaroglu, Ö. (2020) <doi:10.1142/S1088424619500196> Synthesis, characterization and VOCs adsorption kinetics of diethylstilbestrol-substituted metallophthalocyanines. In Porphyrin Science By Women (In 3 Volumes) (pp. 991-999). World Scientific Publishing Co.
#' @export
richie <- function(t,qt,qe,n){
  x <- t
  y <- qt
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(x,y){
    fxn  <- y ~ (qe* (1- exp(-a*x)))
    grd1 <- data.frame(qe = c(0,200),
                        a = c(0,10))
    cc<- capture.output(type="message",
                        fit7 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                     silent=TRUE))
    if(is.null(fit7)==TRUE){
      fit7 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsre <- as.vector(coefficients(fit7))
      pars_qe <- parsre[1L]; pars_a <- parsre[2L]; pars_lin <- parsre[3L]
      qemin = pars_qe*0.9*pars_lin   ;qemax = pars_qe*1.1*pars_lin
      amin = pars_a*0.9              ;amax = pars_a*1.1

      grd2 <- data.frame(qe=c(qemin,qemax),
                         a=c(amin,amax))
      fit7 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val<- predict(fit7)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit7))
      mae  <- mae(y,predict(fit7))
      mse  <- mse(y,predict(fit7))
      rae  <- rae(y,predict(fit7))
      PAIC <- AIC(fit7)
      PBIC <- BIC(fit7)
      SE   <- sqrt((sum((y-predict(fit7))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Richie Model n=1")
    print(summary(fit7))
    predval(x,n.dat)
    error(y)
    parsre1 <- as.vector(coefficients(fit7))
    pars_qe <- parsre1[1L]; pars_a <- parsre1[2L]
    theme_set(theme_bw())
    fun.1 <- function(x) (pars_qe* (1- exp(-pars_a*x)))
    plot  <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Richie model n=1",
           y="qt",
           x="time",
           title="Richie Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(x,y){
    fxn <- y ~ ((a*qe*x)/(1 + (a*x)))
    grd1 <- data.frame(qe = c(0,200),
                       a  = c(0,10))
    cc<- capture.output(type="message",
                        fit7 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit7)==TRUE){
      fit7 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control= list(maxiter = 1000))
      parsre  <- as.vector(coefficients(fit7))
      pars_qe <- parsre[1L]; pars_a <- parsre[2L]; pars_lin <- parsre[3L]
      fun.1 <- function(x){((pars_lin*pars_a*pars_qe*x)/(1 + (pars_a*x)))}
      r <- fun.1(100000)
      qemin = r*0.9      ; qemax = r*1.1
      amin = pars_a*0.9  ; amax = pars_a*1.1
      grd2 <- data.frame(qe=c(qemin,qemax),
                         a=c(amin,amax))
      fit7 <- nls2(fxn,
                     start = grd2,
                     algorithm = "brute-force",
                     control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit7)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit7))
      mae  <- mae(y,predict(fit7))
      mse  <- mse(y,predict(fit7))
      rae  <- rae(y,predict(fit7))
      PAIC <- AIC(fit7)
      PBIC <- BIC(fit7)
      SE   <- sqrt((sum((y-predict(fit7))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Richie Model n=2")
    print(summary(fit7))
    predval(x,n.dat)
    error(y)
    parsre1 <- as.vector(coefficients(fit7))
    pars_qe <- parsre1[1L]; pars_a <- parsre1[2L]
    theme_set(theme_bw())
    fun.1 <- function(x) ((pars_a*pars_qe*x)/(1 + (pars_a*x)))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Richie model n=2",
           y="qt",
           x="time",
           title="Richie Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQn <- function(x,y){
    fxn <- y ~ qe-(qe*(1+((n-1)*a*x))^(1/1-n))
    grd1 <- data.frame(qe = c(0,100),
                       a = c(0,100),
                       n = c(1,10))
    cc<- capture.output(type="message",
                        fit7 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter=1000)),
                                    silent=TRUE))
    if(is.null(fit7)==TRUE){
      fit7 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 100))
      parsre <- as.vector(coefficients(fit7))
      pars_qe <- parsre[1L]; pars_a <- parsre[2L];pars_n <- parsre[3L]; pars_lin <- parsre[4L]
      fun.1 <- function(x){(pars_qe*pars_lin - (pars_qe*pars_lin*(1 + (pars_n-1)*pars_a*x))^(1/(1-pars_n)))}
      r <- fun.1(100000000)
      qemin = r*0.9       ;qemax = r*1.1
      amin = pars_a*0.9   ;amax = pars_a*1.1
      nmin = pars_n*0.9   ;nmax = pars_n*1.1
      grd2 <- data.frame(qe=c(qemin,qemax),
                         a=c(amin,amax),
                         n=c(nmin,nmax))
      fit7 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit7)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit7))
      mae  <- mae(y,predict(fit7))
      mse  <- mse(y,predict(fit7))
      rae  <- rae(y,predict(fit7))
      PAIC <- AIC(fit7)
      PBIC <- BIC(fit7)
      SE   <- sqrt((sum((y-predict(fit7))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Richie Model n\U2260\U0031")
    print(summary(fit7))
    predval(x,n.dat)
    error(y)
    parsre1 <- as.vector(coefficients(fit7))
    pars_qe <- parsre1[1L]; pars_a <- parsre1[2L]; pars_n<- parsre1[3L]
    theme_set(theme_bw())
    fun.1 <- function(x) pars_qe-(pars_qe*(1+((pars_n-1)*pars_a*x))^(1/1-pars_n))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Richie model n\U2260\U0031",
           y="qt",
           x="time",
           title="Richie Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQn2 <- function(x,y,n){
    n=n
    fxn <- y ~ qe-(qe*(1+((n-1)*a*x))^(1/1-n))
    grd1 <- data.frame(qe = c(0,100),
                       a = c(0,100))
    cc<- capture.output(type="message",
                        fit7 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter=1000)),
                                    silent=TRUE))
    if(is.null(fit7)==TRUE){
      fit7 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 100))
      parsre <- as.vector(coefficients(fit7))
      pars_qe <- parsre[1L]; pars_a <- parsre[2L]; pars_lin <- parsre[3L]
      fun.1 <- function(x){(pars_qe*pars_lin - (pars_qe*pars_lin*(1 + (n-1)*pars_a*x))^(1/(1-n)))}
      r <- fun.1(100000000)
      qemin = r*0.9;qemax = r*1.1
      amin = pars_a*0.9;amax = pars_a*1.1
      grd2 <- data.frame(qe=c(qemin,qemax),
                           a=c(amin,amax))
      fit7 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit7)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(x){
      rmse <- rmse(x,predict(fit7))
      mae  <- mae(x,predict(fit7))
      mse  <- mse(x,predict(fit7))
      rae  <- rae(x,predict(fit7))
      PAIC <- AIC(fit7)
      PBIC <- BIC(fit7)
      SE   <- sqrt((sum((y-predict(fit7))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Richie Model n=",n,sep="")
    print(summary(fit7))
    predval(x,n.dat)
    error(y)
    numb<- n
    parsre1 <- as.vector(coefficients(fit7))
    pars_qe <- parsre1[1L]; pars_a <- parsre1[2L]
    theme_set(theme_bw())  # pre-set the bw theme.
    fun.1 <- function(x) pars_qe-(pars_qe*(1+((n-1)*pars_a*x))^(1/1-n))
    subtitle_input <- capture.output(message("Plot of time vs qt with non-linear Richie model n=",numb,sep=""),type="message")
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle=subtitle_input,
           y="qt",
           x="time",
           title="Richie Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ1qe <- function(x,y,qe){
    qe <- qe
    fxn  <- y ~ (qe* (1- exp(-a*x)))
    grd1 <- data.frame(a = c(0,10))
    cc<- capture.output(type="message",
                        fit7 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit7)==TRUE){
      fit7 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsre <- as.vector(coefficients(fit7))
      pars_a <- parsre[1L]
      amin = pars_a*0.9   ;amax = pars_a*1.1
      grd2 <- data.frame(a=c(amin,amax))
      fit7 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val<- predict(fit7)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit7))
      mae  <- mae(y,predict(fit7))
      mse  <- mse(y,predict(fit7))
      rae  <- rae(y,predict(fit7))
      PAIC <- AIC(fit7)
      PBIC <- BIC(fit7)
      SE   <- sqrt((sum((y-predict(fit7))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Richie Model n=1")
    print(summary(fit7))
    predval(x,n.dat)
    error(y)
    parsre1 <- as.vector(coefficients(fit7))
    pars_a <- parsre1[1L]
    theme_set(theme_bw())
    fun.1 <- function(x) (qe* (1- exp(-pars_a*x)))
    plot  <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of time vs qt with non-linear Richie model n=1",
           y="qt",
           x="time",
           title="Richie Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2qe <- function(x,y,qe){
    qe  <- qe
    fxn <- y ~ ((a*qe*x)/(1 + (a*x)))
    grd1 <- data.frame(a  = c(0,10))
    cc<- capture.output(type="message",
                        fit7 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit7)==TRUE){
      fit7 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control= list(maxiter = 1000))
      parsre  <- as.vector(coefficients(fit7))
      pars_a <- parsre[1L]
      amin = pars_a*0.9  ; amax = pars_a*1.1
      grd2 <- data.frame(a=c(amin,amax))
      fit7 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit7)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit7))
      mae  <- mae(y,predict(fit7))
      mse  <- mse(y,predict(fit7))
      rae  <- rae(y,predict(fit7))
      PAIC <- AIC(fit7)
      PBIC <- BIC(fit7)
      SE   <- sqrt((sum((y-predict(fit7))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Richie Model n=2")
    print(summary(fit7))
    predval(x,n.dat)
    error(y)
    parsre1 <- as.vector(coefficients(fit7))
    pars_a <- parsre1[1L]
    theme_set(theme_bw())
    fun.1 <- function(x) ((pars_a*qe*x)/(1 + (pars_a*x)))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of time vs qt with non-linear Richie model n=2",
           y="qt",
           x="time",
           title="Richie Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQnqe <- function(x,y,qe){
    qe  <- qe
    fxn <- y ~ qe-(qe*(1+((n-1)*a*x))^(1/1-n))
    grd1 <- data.frame(a = c(0,100),
                       n = c(1,10))
    cc<- capture.output(type="message",
                        fit7 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter=1000)),
                                    silent=TRUE))
    if(is.null(fit7)==TRUE){
      fit7 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 100))
      parsre <- as.vector(coefficients(fit7))
      pars_a <- parsre[1L];pars_n <- parsre[2L]
      amin = pars_a*0.9   ;amax = pars_a*1.1
      nmin = pars_n*0.9   ;nmax = pars_n*1.1
      grd2 <- data.frame(a=c(amin,amax),
                         n=c(nmin,nmax))
      fit7 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit7)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit7))
      mae  <- mae(y,predict(fit7))
      mse  <- mse(y,predict(fit7))
      rae  <- rae(y,predict(fit7))
      PAIC <- AIC(fit7)
      PBIC <- BIC(fit7)
      SE   <- sqrt((sum((y-predict(fit7))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Richie Model n\U2260\U0031")
    print(summary(fit7))
    predval(x,n.dat)
    error(y)
    parsre1 <- as.vector(coefficients(fit7))
    pars_a <- parsre1[1L]; pars_n<- parsre1[2L]
    theme_set(theme_bw())
    fun.1 <- function(x) qe-(qe*(1+((pars_n-1)*pars_a*x))^(1/1-pars_n))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of time vs qt with non-linear Richie model n\U2260\U0031",
           y="qt",
           x="time",
           title="Richie Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQn2qe <-function(x,y,qe,n){
    n   <- n
    qe  <- qe
    fxn <- y ~ qe-(qe*(1+((n-1)*a*x))^(1/1-n))
    grd1 <- data.frame(a = c(0,100))
    cc<- capture.output(type="message",
                        fit7 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter=1000)),
                                    silent=TRUE))
    if(is.null(fit7)==TRUE){
      fit7 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 100))
      parsre <- as.vector(coefficients(fit7))
      pars_a <- parsre[1L]
      amin = pars_a*0.9; amax = pars_a*1.1
      grd2 <- data.frame(a=c(amin,amax))
      fit7 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit7)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(x){
      rmse <- rmse(x,predict(fit7))
      mae  <- mae(x,predict(fit7))
      mse  <- mse(x,predict(fit7))
      rae  <- rae(x,predict(fit7))
      PAIC <- AIC(fit7)
      PBIC <- BIC(fit7)
      SE   <- sqrt((sum((y-predict(fit7))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Richie Model n=",n,sep="")
    print(summary(fit7))
    predval(x,n.dat)
    error(y)
    numb<- n
    parsre1 <- as.vector(coefficients(fit7))
    pars_a <- parsre1[1L]
    theme_set(theme_bw())
    fun.1 <- function(x) qe-(qe*(1+((n-1)*pars_a*x))^(1/1-n))
    subtitle_input <- capture.output(message("Plot of qt vs time with non-linear Richie model n=",numb,sep=""),type="message")
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle=subtitle_input,
           y="qt",
           x="time",
           title="Richie Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  if(missing(qe)){
    if(missing(n)){
      EQn(x,y)}
    else if(is.null(n)){
      EQn(x,y)}
    else if(isFALSE(n)){
      EQn(x,y)}
    else if(n==1){
      EQ1(x,y)}
    else if(n==2){
      EQ2(x,y)}
    else{
      EQn2(x,y,n)}
    }
  else if(is.null(qe)){
    if(missing(n)){
      EQn(x,y)}
    else if(is.null(n)){
      EQn(x,y)}
    else if(isFALSE(n)){
      EQn(x,y)}
    else if(n==1){
      EQ1(x,y)}
    else if(n==2){
      EQ2(x,y)}
    else{
      EQn2(x,y,n)}}
  else if(isFALSE(qe)){
    if(missing(n)){
      EQn(x,y)}
    else if(is.null(n)){
      EQn(x,y)}
    else if(isFALSE(n)){
      EQn(x,y)}
    else if(n==1){
      EQ1(x,y)}
    else if(n==2){
      EQ2(x,y)}
    else{
      EQn2(x,y,n)}}
  else{
    if(missing(n)){
      EQnqe(x,y,qe)}
    else if(is.null(n)){
      EQnqe(x,y,qe)}
    else if(isFALSE(n)){
      EQnqe(x,y,qe)}
    else if(n==1){
      EQ1qe(x,y,qe)}
    else if(n==2){
      EQ2qe(x,y,qe)}
    else{
      EQn2qe(x,y,qe,n)}
  }
}
