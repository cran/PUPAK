#' @title Non-Linear Pseudo-First-Order Adsorption Kinetic Model
#' @description The Pseudo-First Order Adsorption Kinetic Model follows the Linear Driving Force model (LDF) which states that the rate of mass transfer is equal to the transfer coefficient and the difference between the amount adsorbed and the amount adsorbed in equilibrium. This model is an empirical rate equation known to describe the rate of sorption in liquid-phase systems. The PFO model is not suitable for the whole adsorption reaction since the rate of adsorption decreases until it reaches the maximum adsorption capacity, and thus, the rate is zero at equilibrium (Plazinski, Rudzinski, and Plazinska, 2009).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qe the numerical value for the amount adsorbed at equilibrium. If this parameter is not defined, it will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @import utils
#' @return the non-linear regression and the parameter estimation for the Pseudo-First-Order Model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples qe <- 4.68
#' @examples pfo(t,qt,qe)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Lagergren, S. (1898), Zur theorie der sogenannten adsorption gelöster stoffe, Kungliga Svenska Vetenskapsakademiens. Handlingar, 24 (4) : 1-39.
#' @references Plazinski, W., Rudzinski, W., &; Plazinska, A. (2009) <doi:10.1016/j.cis.2009.07.009> Theoretical models of sorption kinetics including a surface reaction mechanism: A review. In Advances in Colloid and Interface Science (Vol. 152, Issues 1-2, pp. 2-13).
#' @export
pfo <- function(t,qt,qe){
  x   <- t
  y   <- qt
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(x,y){
    fxn <- y ~ (qe* (1- exp(-k1*x)))
    grd1 <- data.frame(qe = c(0,1000),
                       k1 = c(0,10))
    cc<- capture.output(type="message",
                        fit1 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.character(fit1)==TRUE){
      fit1 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parspfo <- as.vector(coefficients(fit1))
      pars_qe <- parspfo[1L]; pars_k1 <- parspfo[2L]; pars_lin <- parspfo[3L]
      qemin <- pars_qe*0.9*pars_lin ; qemax <- pars_qe*1.1*pars_lin
      k1min <- pars_k1*0.9          ; k1max <- pars_k1*1.1
      grd2  <- data.frame(qe=c(qemin,qemax),
                          k1=c(k1min,k1max))
      fit1 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit1)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit1))
      mae  <- mae(y,predict(fit1))
      mse  <- mse(y,predict(fit1))
      rae  <- rae(y,predict(fit1))
      PAIC <- AIC(fit1)
      PBIC <- BIC(fit1)
      SE   <- sqrt((sum((y-predict(fit1))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Pseudo First Order Model")
    print(summary(fit1))
    predval(x,n.dat)
    error(y)
    parspfo1 <- as.vector(coefficients(fit1))
    pars_qe <- parspfo1[1L]; pars_k1 <- parspfo1[2L]
    theme_set(theme_bw())
    fun.1 <- function(x) (pars_qe* (1- exp(-pars_k1*x)))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Pseudo First Order Model",
           y="qt",
           x="time",
           title="Pseudo First Order Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(x,y,qe){
    qe <-qe
    fxn <- y ~ (qe* (1- exp(-k1*x)))
    grd1 <- data.frame(k1 = c(0,10))
    cc<- capture.output(type="message",
                        fit1 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1500)),
                                    silent=TRUE))
    if(is.character(fit1)==TRUE){
      fit1 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parspfo <- as.vector(coefficients(fit1))
      pars_k1 <- parspfo[1L]
      k1min <- pars_k1*0.9; k1max <- pars_k1*1.1
      grd2  <- data.frame(k1=c(k1min,k1max))
      fit1 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit1)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit1))
      mae  <- mae(y,predict(fit1))
      mse  <- mse(y,predict(fit1))
      rae  <- rae(y,predict(fit1))
      PAIC <- AIC(fit1)
      PBIC <- BIC(fit1)
      SE   <- sqrt((sum((y-predict(fit1))^2))/(n.dat-2))

      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Pseudo First Order Model")
    print(summary(fit1))
    predval(x,n.dat)
    error(y)
    parspfo1 <- as.vector(coefficients(fit1))
    pars_k1 <- parspfo1[1L]
    theme_set(theme_bw())
    fun.1 <- function(x) (qe* (1- exp(-pars_k1*x)))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Pseudo First Order Model",
           y="qt",
           x="time",
           title="Pseudo First Order Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  if(missing(qe)){
    EQ1(x,y)}
  else if(is.null(qe)){
    EQ1(x,y)}
  else if(isFALSE(qe)){
    EQ1(x,y)}
  else{
    EQ2(x,y,qe)}
}
