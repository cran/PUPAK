#' @title Non-Linear Elovich Adsorption Kinetic Model
#' @description The Elovich Adsorption Kinetic Model is an empirical rate equation that states the adsorption energy rises in a linear relationship with surface coverage. The model assumes that adsorption occurs on localized sites, the interaction between adsorbed ions is present, and the concentration of adsorbate is considered to be constant. It is applicable in gas adsorptions as well as wastewater processes (Largitte and Pasquier, 2016).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameter estimation for the Elovich Adsorption Kinetic Model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' elovich(t,qt)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Elovich, S. Y., & Larinov, O. G. (1962) Theory of adsorption from solutions of non electrolytes on solid (I) equation adsorption from solutions and the analysis of its simplest form,(II) verification of the equation of adsorption isotherm from solutions. Izvestiya Akademii Nauk, 2(2), 209–216.
#' @references Largitte, L., &; Pasquier, R. (2016) <doi:10.1016/j.cherd.2016.02.006> A review of the kinetics adsorption models and their application to the adsorption of lead by an activated carbon. Chemical Engineering Research and Design, 109, 495-504.
#' @export
elovich <- function(t,qt){
  x <- t
  y <- qt
  fxn <- y ~ ((1/beta)* log(1 + alpha*beta*x))
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  grd1 <- data.frame(alpha = c(0,1000),
                     beta = c(0,10))
  cc<- capture.output(type="message",
                      fit5 <- try(nls2::nls2(fxn,
                                             data = dat,
                                             start = grd1,
                                             algorithm = "port",
                                             control = list(maxiter = 2000)),
                                   silent=TRUE))
  if(is.null(fit5)==TRUE){
    fitem <- nls2(fxn,
                  data = dat,
                  start = grd1,
                  algorithm = "plinear-random",
                  control = list(maxiter = 1000))
    parsem <- as.vector(coefficients(fitem))
    pars_alpha <- parsem[1L]; pars_beta <- parsem[2L]; pars_lin <- parsem[3L]
    alphamin <- pars_alpha*0.9           ;alphamax <- pars_alpha*1.1
    betamin  <- (pars_beta/pars_lin)*0.9 ;betamax  <- (pars_beta/pars_lin)*1.1
    grd2 <- data.frame(alpha=c(alphamin,alphamax),
                         beta=c(betamin,betamax))
    fit5 <- nls2(fxn,
                 start = grd2,
                 algorithm = "brute-force",
                 control=list(maxiter=1000))
  }else{}
  predval <- function(x,n.dat){
    Col1<- c(rep(" |",each = n.dat))
    Col2<- c(rep("|",each = n.dat))
    pred.val <- predict(fit5)
    time <- x
    P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
    colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
    message("Estimated Values")
    print(P.Table, right=T, row.names = F)
  }
  error <- function(y){
    rmse <- (rmse(y,predict(fit5)))
    mae  <- (mae(y,predict(fit5)))
    mse  <- (mse(y,predict(fit5)))
    rae  <- (rae(y,predict(fit5)))
    PAIC <- AIC(fit5)
    PBIC <- BIC(fit5)
    SE   <- sqrt((sum((y-predict(fit5))^2))/(n.dat-2))
    Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
    Col2 <- c("|","|","|","|","|","|","|")
    E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
    E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
    E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
    colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
    message("Error Estimation")
    print(E.Table, right=F, row.names = F)
  }
  message("Elovich Adsorption Kinetic Model")
  print(summary(fit5))
  predval(x,n.dat)
  error(y)
  parsem <- as.vector(coefficients(fit5))
  pars_alpha <- parsem[1L]; pars_beta <- parsem[2L]
  theme_set(theme_bw())
  fun.1 <- function(x){(1/pars_beta)* log(1 + pars_alpha*pars_beta*x)}
  plot <- ggplot(dat, aes(x=x,y=y))+
    geom_function(color="red", fun=fun.1,size=1)+
    geom_point()+
    labs(subtitle="Plot of qt vs time with non-linear Elovich model",
         y="qt",
         x="time",
         title="Elovich Model",
         caption="Created by PUPAK using ggplot2")
  print(plot)
}
