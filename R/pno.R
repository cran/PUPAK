#' @title Non-Linear Pseudo-nth-Order Adsorption Kinetic Model
#' @description The Pseudo-nth Order Adsorption Kinetic Model is an empirical rate equation known to describe the kinetic analysis of neither order 1 nor 2 kinetic parameters. It will have a significant effect on the calculation of the rate constants as the rate constant is dependent on the order of reaction (Tseng, Wu, Wu, and Juang, 2014).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qe the numerical value for the amount adsorbed at equilibrium. If this parameter is not defined, it will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameter estimation for the Pseudo-nth-Order Model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qe <- 4.8
#' pno(t,qt,qe)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Ozer, A. (2007) <doi:10.1016/j.jhazmat.2006.07.040> Removal of Pb(II) ions from aqueous solutions by sulphuric acid-treated wheat bran. Journal of Hazardous Materials, 141(3), 753-761.
#' @references Tseng, R. L., Wu, P. H., Wu, F. C., &; Juang, R. S. (2014) <doi:10.1016/j.cej.2013.10.013> A convenient method to determine kinetic parameters of adsorption processes by nonlinear regression of pseudo-nth-order equation. Chemical Engineering Journal, 237, 153-161.
#' @export
pno <- function(t,qt,qe){
  x <- t
  y <- qt
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(x,y){
    fxn <- y ~ qe*(1-(1/((1+((n-1)*kn*x*(qe^(n-1))))^(1/(n-1)))))
    grd1  <- data.frame(qe = c(0,1000),
                        kn = c(0,10),
                        n = c(1,3))
    cc<- capture.output(type="message",
                        fit3 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "plinear-random",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    pars <- as.vector(coefficients(fit3))
    pars_qe <- pars[1L]; pars_kn <- pars[2L]; pars_n <- pars[3L]; pars_lin <- pars[4L]
    fun.1 <- function(x)(pars_qe*pars_lin*(1-(1/((1+((pars_n-1)*pars_kn*x*(pars_qe^(pars_n-1))))^(1/(pars_n-1))))))
    r <- fun.1(1000)
    qemin <- r*0.9       ;qemax <- r*1.1
    knmin <- pars_kn*0.9 ;knmax <- pars_kn*1.1
    nmin  <- pars_n*0.9  ;nmax  <- pars_n*1.1
    grd2 <- data.frame(qe=c(qemin,qemax),
                       kn=c(knmin,knmax),
                       n =c(nmin,nmax))
    fit3 <- nls2(fxn,
                 start = grd2,
                 algorithm = "grid-search",
                 control=list(maxiter=1000))

    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit3)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit3))
      mae  <- mae(y,predict(fit3))
      mse  <- mse(y,predict(fit3))
      rae  <- rae(y,predict(fit3))
      PAIC <- AIC(fit3)
      PBIC <- BIC(fit3)
      SE   <- sqrt((sum((y-predict(fit3))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Pseudo nth Order")
    print(summary(fit3))
    predval(x,n.dat)
    error(y)
    parspno1 <- as.vector(coefficients(fit3))
    pars_qe <- parspno1[1L]; pars_kn <- parspno1[2L]; pars_n <- parspno1[3L]
    theme_set(theme_bw())
    fun.1 <- function(x) (pars_qe*(1-(1/((1+((pars_n-1)*pars_kn*x*(pars_qe^(pars_n-1))))^(1/(pars_n-1))))))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with nonlinear Pseudo nth Order model",
           y="qt",
           x="time",
           title="Pseudo nth Order Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(x,y,qe){
    qe <-qe
    fxn <- y ~ qe*(1-(1/((1+((n-1)*kn*x*(qe^(n-1))))^(1/(n-1)))))
    grd1  <- data.frame(kn = c(0,10),
                        n = c(1,3))
    cc<- capture.output(type="message",
                        fit3 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.character(fit3)==TRUE){
      fit3 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      pars <- as.vector(coefficients(fit3))
      pars_kn <- pars[1L]; pars_n <- pars[2L]
      knmin <- pars_kn*0.9 ;knmax <- pars_kn*1.1
      nmin  <- pars_n*0.9  ;nmax  <- pars_n*1.1
      grd2 <- data.frame(kn=c(knmin,knmax),
                         n =c(nmin,nmax))
      fit3 <- nls2(fxn,
                   start = grd2,
                   algorithm = "grid-search",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit3)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit3))
      mae  <- mae(y,predict(fit3))
      mse  <- mse(y,predict(fit3))
      rae  <- rae(y,predict(fit3))
      PAIC <- AIC(fit3)
      PBIC <- BIC(fit3)
      SE   <- sqrt((sum((y-predict(fit3))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Pseudo nth Order")
    print(summary(fit3))
    predval(x,n.dat)
    error(y)
    parspno1 <- as.vector(coefficients(fit3))
    pars_kn <- parspno1[1L]; pars_n <- parspno1[2L]
    theme_set(theme_bw())
    fun.1 <- function(x) (qe*(1-(1/((1+((pars_n-1)*pars_kn*x*(qe^(pars_n-1))))^(1/(pars_n-1))))))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with nonlinear Pseudo nth Order model",
           y="qt",
           x="time",
           title="Pseudo nth Order Model",
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
