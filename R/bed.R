#' @title Non-Linear Boyd External Diffusion Model
#' @description The Boyd External Diffusion Model is frequently applied to adsorption kinetic data to calculate the rate constant, assuming that film diffusion is the rate-limiting step in the first few minutes of the adsorption process. The film diffusion has a strong dependency on agitation. Boyd’s diffusion models are used in numerous adsorption studies mostly to determine the rate-controlling step (Viegas, Campinas, Costa, and Rosa, 2014).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qinf the numerical value for the amount adsorbed at infinite time. If this parameter is not defined, it will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameter estimation for the Boyd External Diffusion model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qinf <- 4.68
#' bed(t,qt,qinf)
#' bed(t,qt)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Boyd, G. E., Adamson, A. W., & Myers, L. S. (1947) <doi:10.1021/ja01203a066> The Exchange Adsorption of Ions from Aqueous Solutions by Organic Zeolites. II. Kinetics1. Journal of the American Chemical Society, 69(11), 2836-2848.
#' @references Viegas, R. M. C., Campinas, M., Costa, H., &; Rosa, M. J. (2014) <doi:10.1007/s10450-014-9617-9> How do the HSDM and Boyd's model compare for estimating intraparticle diffusion coefficients in adsorption processes. Adsorption, 20(5-6), 737-746.
#' @export
bed <- function(t,qt,qinf){
  x <- t
  y <- qt
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(t,qt){
    fxn  <- y ~ (qinf*(1 - exp(-R*x)))
    grd1 <- data.frame(qinf = c(0,1000),
                       R    = c(0,10))
    cc<- capture.output(type="message",
                        fit8 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit8)==TRUE){
      fit8 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsbed <- as.vector(coefficients(fit8))
      pars_qinf <- parsbed[1L]; pars_R <- parsbed[2L]; pars_lin <- parsbed[3L]
      qinfmin <- pars_qinf*0.9*pars_lin    ; qinfmax <- pars_qinf*1.1*pars_lin
      Rmin    <- pars_R*0.9                ; Rmax <- pars_R*1.1
      grd2 <- data.frame(qinf=c(qinfmin,qinfmax),
                         R=c(Rmin,Rmax))
      fit8 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit8)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- (rmse(y,predict(fit8)))
      mae  <- (mae(y,predict(fit8)))
      mse  <- (mse(y,predict(fit8)))
      rae  <- (rae(y,predict(fit8)))
      PAIC <- AIC(fit8)
      PBIC <- BIC(fit8)
      SE   <- sqrt((sum((y-predict(fit8))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Boyd External Diffusion Model")
    print(summary(fit8))
    predval(x,n.dat)
    error(y)
    parsbed1 <- as.vector(coefficients(fit8))
    pars_qinf1 <- parsbed1[1L]; pars_R1 <- parsbed1[2L]
    theme_set(theme_bw())
    fun.1 <- function(x) {pars_qinf1*(1 - exp(-pars_R1*x))}
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of time vs qt with non-linear Boyd External Diffusion model",
           y="qt",
           x="time",
           title="Boyd External Diffusion Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(t,qt,qinf){
    qinf <- qinf
    fxn  <- y ~ (qinf*(1 - exp(-R*x)))
    grd1 <- data.frame(R = c(0,10))
    cc<- capture.output(type="message",
                        fit8 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit8)==TRUE){
      fit8 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsbed <- as.vector(coefficients(fit8))
      pars_R <- parsbed[1L]; pars_lin <- parsbed[2L]
      Rmin <- pars_R*0.9  ;Rmax <- pars_R*1.1
      grd2 <- data.frame(R=c(Rmin,Rmax))
      fit8 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit8)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- (rmse(y,predict(fit8)))
      mae  <- (mae(y,predict(fit8)))
      mse  <- (mse(y,predict(fit8)))
      rae  <- (rae(y,predict(fit8)))
      PAIC <- AIC(fit8)
      PBIC <- BIC(fit8)
      SE   <- sqrt((sum((y-predict(fit8))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Boyd External Diffusion Model")
    print(summary(fit8))
    predval(x,n.dat)
    error(y)
    parsbed1 <- as.vector(coefficients(fit8))
    pars_R1  <- parsbed1[1L]
    theme_set(theme_bw())
    fun.1 <- function(x) {qinf*(1 - exp(-pars_R1*x))}
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with nonlinear Boyd External Diffusion model",
           y="qt",
           x="time",
           title="Boyd External Diffusion Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  if(missing(qinf)){
    EQ1(x,y)}
  else if(is.null(qinf)){
    EQ1(x,y)}
  else if(isFALSE(qinf)){
    EQ1(x,y)}
  else{EQ2(x,y,qinf)
    }
}
