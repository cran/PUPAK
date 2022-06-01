#' @title Non-Linear Matthews and Weber Diffusion Model
#' @description The Matthews and Weber Diffusion Model is used in determining the rate constant for film diffusion where it assumes that intraparticle diffusion can be neglected at the early period of contact. Derived from Fickien’s law application, the solute concentration in the liquid phase in this model is expressed as a function of solute concentration difference in the liquid phase and at the adsorbent surface (Prasad & Srivastava, 2009).
#' @param t the numerical value for contact time
#' @param Ct the numerical value for the concentration of the adsorbent at time t
#' @param Co the numerical value for the initial concentration of the adsorbent. If this parameter is not defined, it will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameter estimation for the Matthews and Weber diffusion model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' Ct <- c(10.000,8.141,8.056,7.949,7.863,7.799,7.778,7.756,7.692)
#' Co <- 10
#' maw(t,Ct,Co)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Mathews, A. P., &; Weber, W. J. (1984) <doi:10.1080/00986448408940104> Modeling and parameter evaluation for adsorption in slurry reactors. Chemical Engineering Communications, 25(1-6), 157-171.
#' @references Krishna Prasad, R., & Srivastava, S. N. (2009) <doi:10.1016/j.cej.2008.05.021> Sorption of distillery spent wash onto fly ash: Kinetics and mass transfer studies. Chemical Engineering Journal, 146(1), 90–97.
#' @export

maw <- function(t,Ct,Co){
  x <- t; y <- Ct
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1<-function(x,y){
    fxn <- y ~ Co*exp(-kmaw*x)
    grd1 <- data.frame(Co = c(0,100),
                       kmaw  = c(0,1))
    cc<- capture.output(type="message",
                        fit10 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit10)==TRUE){
      fit10 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsmaw <- as.vector(coefficients(fit10))
      pars_Co <- parsmaw[1L]; pars_kmaw <- parsmaw[2L]; pars_lin <- parsmaw[3L]
      fun.1 <- function(x)(pars_Co*pars_lin*exp(-pars_kmaw*x))
      r <- fun.1(0)
      Comin <- r*0.9            ;Comax   <- r*1.1
      kmawmin <-pars_kmaw*0.9   ;kmawmax <-pars_kmaw*1.1
      grd2 <- data.frame(Co=c(Comin,Comax),
                         kmaw=c(kmawmin,kmawmax))
      fit10 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit10)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit10))
      mae  <- mae(y,predict(fit10))
      mse  <- mse(y,predict(fit10))
      rae  <- rae(y,predict(fit10))
      PAIC <- AIC(fit10)
      PBIC <- BIC(fit10)
      SE   <- sqrt((sum((y-predict(fit10))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Matthews and Weber Diffusion Model")
    print(summary(fit10))
    predval(x,n.dat)
    error(y)

    parsmaw <- as.vector(coefficients(fit10))
    pars_Co <- parsmaw[1L]; pars_kmaw <- parsmaw[2L]
    theme_set(theme_bw())
    fun.1 <- function(x)(pars_Co*exp(-pars_kmaw*x))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Matthews and Weber model",
           y="qt",
           x="time",
           title="Matthews and Weber Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2<-function(x,y,Co){
    fxn <- y ~ Co*exp(-kmaw*x)
    grd1 <- data.frame(kmaw  = c(0,1))
    cc<- capture.output(type="message",
                        fit10 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit10)==TRUE){
      fit10 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsmaw <- as.vector(coefficients(fit10))
      kmawmin =  pars_kmaw*0.9  ;kmawmax =  pars_kmaw*1.1
      grd2 <- data.frame(kmaw=c(kmawmin,kmawmax))
      fit10 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit10)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit10))
      mae  <- mae(y,predict(fit10))
      mse  <- mse(y,predict(fit10))
      rae  <- rae(y,predict(fit10))
      PAIC <- AIC(fit10)
      PBIC <- BIC(fit10)
      SE   <- sqrt((sum((y-predict(fit10))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }

    message("Matthews and Weber Diffusion Model")
    print(summary(fit10))
    predval(x,n.dat)
    error(y)

    parsmaw <- as.vector(coefficients(fit10))
    pars_kmaw <- parsmaw[1L]
    theme_set(theme_bw())
    fun.1 <- function(x)(Co*exp(-pars_kmaw*x))
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear Matthews and Weber model",
           y="qt",
           x="time",
           title="Matthews and Weber Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  if(missing(Co)){
    EQ1(x,y)}
  else if(is.null(Co)){
    EQ1(x,y)}
  else if(isFALSE(Co)){
    EQ1(x,y)}
  else{
    EQ2(x,y,Co)}
}
