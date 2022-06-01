#' @title Non-Linear Fractional Power Adsorption Kinetic Model
#' @description The Fractional Power Adsorption Kinetic Model is an empirical rate equation in which the specific adsorption rate at a unit time can be estimated using its product's constant.(Netzahuatl-Munoz, Del Carmen Cristiani-Urbina, and Cristiani-Urbina, 2015)
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qe the numerical value for the amount adsorbed at equilibrium. If this parameter is not defined, it will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameter estimation for the Fractional Power Adsorption Kinetic Model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qe <- 4.68
#' fp(t,qt,qe)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Dalai, R. C. (1974) <doi:10.1080/00103627409366531> Desorption of soil phosphate by anion-exchange resin. Communications in Soil Science and Plant Analysis, 5(6), 531-538.
#' @references Netzahuatl-Munoz, A. R., del Carmen Cristiani-Urbina, M., &; Cristiani-Urbina, E. (2015) <doi:10.1371/journal.pone.0137086> Chromium biosorption from Cr(VI) aqueous solutions by Cupressus lusitanica bark: Kinetics, equilibrium and thermodynamic studies. PLoS ONE, 10(9).
#' @export
fp <- function(t,qt,qe){
  x   <- t
  y   <- qt
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(x,y){
    fxn  <- y ~ qe*(alpha*(x)^beta)
    grd1 <- data.frame(qe = c(0,1000),
                       alpha = c(0,1000),
                       beta = c(0,10))
    cc<- capture.output(type="message",
                        fit6 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit6)==TRUE){
      fit6 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsfpm  <- as.vector(coefficients(fit6))
      pars_qe  <- parsfpm[1L]; pars_alpha <- parsfpm[2L]; pars_beta<- parsfpm[3L]; pars_lin <- parsfpm[4L]
      fun.1    <- function(x){pars_qe*pars_lin*(pars_alpha*(x)^pars_beta)}
      r        <- fun.1(100000000)
      qemin    <- r*1.25                             ;qemax    <- r*0.75
      alphamin <- pars_qe*pars_alpha*0.9*pars_lin/r  ;alphamax <- pars_qe*pars_alpha*1.1*pars_lin/r
      betamin  <- pars_beta*0.9                      ;betamax  <- pars_beta*1.1
      grd2 <- data.frame(qe=c(qemin,qemax),
                         alpha=c(alphamin,alphamax),
                         beta=c(betamin,betamax))
      fit6 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1     <- c(rep(" |",each = n.dat))
      Col2     <- c(rep("|",each = n.dat))
      pred.val <- predict(fit6)
      time     <- x
      P.Table  <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <-(rmse(y,predict(fit6)))
      mae  <- (mae(y,predict(fit6)))
      mse  <- (mse(y,predict(fit6)))
      rae  <- (rae(y,predict(fit6)))
      PAIC <- AIC(fit6)
      PBIC <- BIC(fit6)
      SE   <- sqrt((sum((y-predict(fit6))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Fractional Power Model")
    print(summary(fit6))
    predval(x,n.dat)
    error(y)
    parsfpm1 <- as.vector(coefficients(fit6))
    pars_qe  <- parsfpm1[1L]; pars_alpha <- parsfpm1[2L]; pars_beta<- parsfpm1[3L]
    theme_set(theme_bw())
    fun.2 <- function(x) {pars_qe*(pars_alpha*(x)^pars_beta)}
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.2,size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with nonlinear Fractional Power Model",
           y="qt",
           x="time",
           title="Fractional Power Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(x,y,qe){
    qe   <- qe
    fxn  <- y ~ qe*(alpha*(x)^beta)
    grd1 <- data.frame(alpha = c(0,1000),
                       beta = c(0,10))
    cc<- capture.output(type="message",
                        fit6 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit6)==TRUE){
      fit6 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsfpm <- as.vector(coefficients(fit6))
      pars_alpha <- parsfpm[1L]; pars_beta<- parsfpm[2L]; pars_lin <- parsfpm[3L]
      fun.1 <- function(x){qe*pars_lin*(pars_alpha*(x)^pars_beta)}
      r <- fun.1(100000000)
      alphamin <- qe*pars_alpha*0.9*pars_lin/r       ;alphamax <- qe*pars_alpha*1.1*pars_lin/r
      betamin  <- pars_beta*0.9                      ;betamax  <- pars_beta*1.1
      grd2 <- data.frame(alpha=c(alphamin,alphamax),
                         beta=c(betamin,betamax))
      fit6 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1     <- c(rep(" |",each = n.dat))
      Col2     <- c(rep("|",each = n.dat))
      pred.val <- predict(fit6)
      time     <- x
      P.Table  <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- (rmse(y,predict(fit6)))
      mae  <- (mae(y,predict(fit6)))
      mse  <- (mse(y,predict(fit6)))
      rae  <- (rae(y,predict(fit6)))
      PAIC <- AIC(fit6)
      PBIC <- BIC(fit6)
      SE   <- sqrt((sum((y-predict(fit6))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Fractional Power Model")
    print(summary(fit6))
    predval(x,n.dat)
    error(y)
    parsfpm1 <- as.vector(coefficients(fit6))
    pars_alpha <- parsfpm1[1L]; pars_beta<- parsfpm1[2L]
    theme_set(theme_bw())
    fun.2 <- function(x) {qe*(pars_alpha*(x)^pars_beta)}
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.2,size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs t with nonlinear Fractional Power Model",
           y="qt",
           x="time",
           title="Fractional Power Model",
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
