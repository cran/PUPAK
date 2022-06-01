#' @title Non-Linear Crank Diffusion Model
#' @description The Crank Diffusion Model is an equation for homogeneous adsorbate diffusion in a sphere-shaped adsorbent with constant surface diffusivity throughout the particle. It's an exact solution for the "infinite bath" case, in which the sphere starts out empty and the solute concentration at the surface remains constant. Due to the constant surface concentration, external film resistance may be ignored (Qiu, Lv, Pan, Zhang, Zhang, and Zhang, 2009).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qinf the numerical value for the amount adsorbed at infinite time. If this argument is not defined, it will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameter estimation for the Crank Adsorption Kinetic Model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qinf <- 4.68
#' crank(t,qt,qinf)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Crank, J. (1979) <ISBN, 0198534116, 9780198534112>The mathematics of diffusion. Oxford university press.
#' @references Qiu, H., Lv, L., Pan, B. C., Zhang, Q. J., Zhang, W. M., &; Zhang, Q. X. (2009) <doi:10.1631/jzus.A0820524> Critical review in adsorption kinetic models. In Journal of Zhejiang University: Science A (Vol. 10, Issue 5, pp. 716-724).
#' @export
crank <- function(t,qt,qinf){
  x <- t
  y <- qt
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(t,qt){
    fxn <- y~(qinf*(1-((6/(pi^2))*((exp(-Dc*(pi^2)*x/(r^2)))+((exp(-Dc*4*(pi^2)*x/(r^2)))/4)+((exp(-Dc*9*(pi^2)*x/(r^2)))/9)
                                +((exp(-Dc*16*(pi^2)*x/(r^2)))/16)+((exp(-Dc*25*(pi^2)*x/(r^2)))/25)+((exp(-Dc*36*(pi^2)*x/(r^2)))/36)
                                +((exp(-Dc*49*(pi^2)*x/(r^2)))/49)+((exp(-Dc*64*(pi^2)*x/(r^2)))/64)+((exp(-Dc*81*(pi^2)*x/(r^2)))/81)
                                +((exp(-Dc*100*(pi^2)*x/(r^2)))/100)))))
    grd1 <- data.frame(qinf = c(0,1000),
                       Dc = c(0,10),
                       r=c(0,10))
    fit11 <- nls2(fxn,
                 data = dat,
                 start = grd1,
                 algorithm = "plinear-random",
                 control = list(maxiter = 2000))
    parscrnk <- as.vector(coefficients(fit11))
    pars_qinf <- parscrnk[1L];pars_Dc <- parscrnk[2L];pars_r <- parscrnk[3L]; pars_lin <- parscrnk[4L]
    qinfmin <- pars_qinf*0.9*pars_lin  ;qinfmax <- pars_qinf*1.1*pars_lin
    Dcmin <- pars_Dc*0.9               ;Dcmax <- pars_Dc*1.1
    rmin <- pars_r*0.9                 ;rmax <- pars_r*1.1
    grd2 <- data.frame(qinf=c(qinfmin,qinfmax),
                       Dc=c(Dcmin,Dcmax),
                       r=c(rmin,rmax))
    fit11 <- nls2(fxn,
                 start = grd2,
                 algorithm = "brute-force",
                 control=list(maxiter=1000))
    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.val <- predict(fit11)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit11))
      mae  <- mae(y,predict(fit11))
      mse  <- mse(y,predict(fit11))
      rae  <- rae(y,predict(fit11))
      PAIC <- AIC(fit11)
      PBIC <- BIC(fit11)
      SE   <- sqrt((sum((y-predict(fit11))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Crank Model")
    print(summary(fit11))
    predval(x,n.dat)
    error(y)
    parscrnk1 <- as.vector(coefficients(fit11))
    pars_qinf <- parscrnk1[1L];pars_Dc <- parscrnk1[2L];pars_r <- parscrnk1[3L]
    theme_set(theme_bw())
    fun.1 <- function(x){pars_qinf*(1-((6/(pi^2))*((exp(-pars_Dc*(pi^2)*x/(pars_r^2)))+((exp(-pars_Dc*4*(pi^2)*x/(pars_r^2)))/4)+((exp(-pars_Dc*9*(pi^2)*x/(pars_r^2)))/9)
                                              +((exp(-pars_Dc*16*(pi^2)*x/(pars_r^2)))/16)+((exp(-pars_Dc*25*(pi^2)*x/(pars_r^2)))/25)+((exp(-pars_Dc*36*(pi^2)*x/(pars_r^2)))/36)
                                              +((exp(-pars_Dc*49*(pi^2)*x/(pars_r^2)))/49)+((exp(-pars_Dc*64*(pi^2)*x/(pars_r^2)))/64)+((exp(-pars_Dc*81*(pi^2)*x/(pars_r^2)))/81)
                                              +((exp(-pars_Dc*100*(pi^2)*x/(pars_r^2)))/100))))}
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with non-linear crank model",
           y="qt",
           x="time",
           title="Crank Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(t,qt,qinf){
    qinf<- qinf
    fxn <- y~(qinf*(1-((6/(pi^2))*((exp(-Dc*(pi^2)*x/(r^2)))+((exp(-Dc*4*(pi^2)*x/(r^2)))/4)+((exp(-Dc*9*(pi^2)*x/(r^2)))/9)
                                   +((exp(-Dc*16*(pi^2)*x/(r^2)))/16)+((exp(-Dc*25*(pi^2)*x/(r^2)))/25)+((exp(-Dc*36*(pi^2)*x/(r^2)))/36)
                                   +((exp(-Dc*49*(pi^2)*x/(r^2)))/49)+((exp(-Dc*64*(pi^2)*x/(r^2)))/64)+((exp(-Dc*81*(pi^2)*x/(r^2)))/81)
                                   +((exp(-Dc*100*(pi^2)*x/(r^2)))/100)))))
    grd1 <- data.frame(Dc = c(0,10),
                       r=c(0,10))
    fit11 <- nls2(fxn,
                  data = dat,
                  start = grd1,
                  algorithm = "plinear-random",
                  control = list(maxiter = 2000))
    parscrnk <- as.vector(coefficients(fit11))
    pars_Dc <- parscrnk[1L]; pars_r <- parscrnk[2L]; pars_lin <- parscrnk[3L]
    Dcmin <- pars_Dc*0.9  ; Dcmax <- pars_Dc*1.1
    rmin  <- pars_r*0.9   ; rmax  <- pars_r*1.1
    grd2 <- data.frame(Dc=c(Dcmin,Dcmax),
                       r=c(rmin,rmax))
    fit11 <- nls2(fxn,
                  start = grd2,
                  algorithm = "brute-force",
                  control=list(maxiter=1000))

    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit11)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit11))
      mae  <- mae(y,predict(fit11))
      mse  <- mse(y,predict(fit11))
      rae  <- rae(y,predict(fit11))
      PAIC <- AIC(fit11)
      PBIC <- BIC(fit11)
      SE   <- sqrt((sum((y-predict(fit11))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Crank Model")
    print(summary(fit11))
    predval(x,n.dat)
    error(y)
    parscrnk1 <- as.vector(coefficients(fit11))
    pars_Dc <- parscrnk1[1L]; pars_r <- parscrnk1[2L]
    theme_set(theme_bw())
    fun.1 <- function(x){qinf*(1-((6/(pi^2))*((exp(-pars_Dc*(pi^2)*x/(pars_r^2)))+((exp(-pars_Dc*4*(pi^2)*x/(pars_r^2)))/4)+((exp(-pars_Dc*9*(pi^2)*x/(pars_r^2)))/9)
                                              +((exp(-pars_Dc*16*(pi^2)*x/(pars_r^2)))/16)+((exp(-pars_Dc*25*(pi^2)*x/(pars_r^2)))/25)+((exp(-pars_Dc*36*(pi^2)*x/(pars_r^2)))/36)
                                              +((exp(-pars_Dc*49*(pi^2)*x/(pars_r^2)))/49)+((exp(-pars_Dc*64*(pi^2)*x/(pars_r^2)))/64)+((exp(-pars_Dc*81*(pi^2)*x/(pars_r^2)))/81)
                                              +((exp(-pars_Dc*100*(pi^2)*x/(pars_r^2)))/100))))}
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1, size=1)+
      geom_point()+
      labs(subtitle="Plot of qt vs time with nonlinear crank model",
           y="qt",
           x="time",
           title="Crank Model",
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
