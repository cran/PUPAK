#' @title Non-Linear Pseudo Second-Order Adsorption Kinetic Model
#' @description The Pseudo-Second-Order Adsorption Kinetic Model is an empirical rate equation known to be the simplified second-order expression of the Pseudo-First Order Adsorption Kinetic Model. It is widely applied to adsorption systems, from biomass to nanomaterials as adsorbent and from heavy metals to pharmaceuticals as adsorbate or contaminant (Revellame, Fortela, Sharp, Hernandez, and Zappi, 2020).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qe the numerical value for the amount adsorbed at equilibrium. If this parameter is not defined, it will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non linear regression and the parameters for the pseudo-second order non-linear model analysis
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qe <- 4.68
#' pso(t,qt,qe)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Ho, Y. S., &; Mckay, G. (1999) <doi:10.1016/S0032-9592(98)00112-5> Pseudo-second order model for sorption processes. In Process Biochemistry (Vol. 34).
#' @references Revellame, E. D., Fortela, D. L., Sharp, W., Hernandez, R., &; Zappi, M. E. (2020) <doi:10.1016/j.clet.2020.100032> Adsorption kinetic modeling using pseudo-first order and pseudo-second order rate laws: A review. In Cleaner Engineering and Technology (Vol. 1). Elsevier Ltd.
#' @export
pso <- function(t,qt,qe){
  x   <- t
  y   <- qt
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(x,y){
    fxn <- y ~ (((qe^2)*k2*x)/(1+(qe*k2*x)))
    grd1 <- data.frame(qe = c(0,1000),
                       k2 = c(0,10))
    cc<- capture.output(type="message",
                        fit2 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.character(fit2)==TRUE){
      fit2 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parspso <- as.vector(coefficients(fit2))
      pars_qe <- parspso[1L]; pars_k2<- parspso[2L]; pars_lin <- parspso[3L]
      fun.1 <- function(x)(((pars_qe^2)*pars_k2*x*pars_lin)/(1+(pars_qe*pars_k2*x)))
      r <- fun.1(100000)
      qemin = r*0.9       ; qemax = r*1.1
      k2min = pars_k2*0.9 ; k2max = pars_k2*1.1
      grd1 <- data.frame(qe=c(qemin,qemax),
                         k2=c(k2min,k2max))
      fit2 <- nls2(fxn,
                   start = grd1,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit2)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit2))
      mae  <- mae(y,predict(fit2))
      mse  <- mse(y,predict(fit2))
      rae  <- rae(y,predict(fit2))
      PAIC <- AIC(fit2)
      PBIC <- BIC(fit2)
      SE   <- sqrt((sum((y-predict(fit2))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Pseudo Second Order")
    print(summary(fit2))
    predval(x,n.dat)
    error(y)
    parspso1 <- as.vector(coefficients(fit2))
    pars_qe  <- parspso1[1L]; pars_k2 <- parspso1[2L]
    theme_set(theme_bw())
    fun.1 <- function(x) (((pars_qe^2)*pars_k2*x)/(1+(pars_qe*pars_k2*x)))
    plot  <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Pseudo Second Order model",
           y="qt",
           x="time",
           title="Pseudo Second Order Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(x,y,qe){
    fxn <- y ~ (((qe^2)*k2*x)/(1+(qe*k2*x)))
    grd1 <- data.frame(k2 = c(0,10))
    cc<- capture.output(type="message",
                        fit2 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.character(fit2)==TRUE){
      fit2 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parspso <- as.vector(coefficients(fit2))
      pars_k2<- parspso[1L]
      k2min = pars_k2*0.9 ; k2max = pars_k2*1.1
      grd1 <- data.frame(k2=c(k2min,k2max))
      fit2 <- nls2(fxn,
                   start = grd1,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit2)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit2))
      mae  <- mae(y,predict(fit2))
      mse  <- mse(y,predict(fit2))
      rae  <- rae(y,predict(fit2))
      PAIC <- AIC(fit2)
      PBIC <- BIC(fit2)
      SE   <- sqrt((sum((y-predict(fit2))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Pseudo Second Order")
    print(summary(fit2))
    predval(x,n.dat)
    error(y)
    parspso1 <- as.vector(coefficients(fit2))
    pars_k2 <- parspso1[1L]
    theme_set(theme_bw())
    fun.1 <- function(x) (((qe^2)*pars_k2*x)/(1+(qe*pars_k2*x)))
    plot  <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Pseudo Second Order model",
           y="qt",
           x="time",
           title="Pseudo Second Order Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  if(missing(qe)){
    EQ1(x,y)}
  else if(is.null(qe)){
    EQ1(x,y)}
  else if(isFALSE(qe)){
    EQ1(x,y)}
  else{EQ2(x,y,qe)}
}
