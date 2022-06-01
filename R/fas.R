#' @title Non-Linear Furusawa and Smith Diffusion Model
#' @description The Furusawa and Smith Diffusion Model is known to describe the rate of adsorption assuming that only external diffusion resistance was predominant during the initial sorption period and controlled the sorption rate. The diffusion model relates the change in fluid phase concentration and time with the fluid phase concentration at the external surface and an external mass transfer coefficient (Furusawa & Smith, 1974).
#' @param t the numerical value for contact time
#' @param Ct the numerical value for the concentration of the adsorbent at time t
#' @param m the numerical value for mass of adsorbent
#' @param V the numerical value for volume of solution
#' @param b the numerical value for the Langmuir isotherm constant
#' @param Co the numerical value for the initial concentration of the adsorbent. If this parameter is not defined, it will be estimated.
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameter estimation for the Furusawa and Smith Diffusion Model
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' Ct <- c(10.000,8.141,8.056,7.949,7.863,7.799,7.778,7.756,7.692)
#' m <- 0.05
#' V <- 0.1
#' b <- 1.3
#' Co <- 10
#' fas(t,Ct,m,V,b,Co)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Furusawa, T., & Smith, J. M. (1974) <doi:10.1002/aic.690200111> Intraparticle mass transport in slurries by dynamic adsorption studies. AIChE Journal, 20(1), 88–93.
#' @references Yakub, E., Agarry, S. E., Omoruwou, F., &; Owabor, C. N. (2020) <doi:10.1080/02726351.2019.1616862> Comparative study of the batch adsorption kinetics and mass transfer in phenol-sand and phenol-clay adsorption systems. Particulate Science and Technology, 38(7), 801-811.
#' @export
fas <- function(t,Ct,m,V,b,Co){
  x <- t     ;y <- Ct
  m <- m     ;V <- V
  dat <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  EQ1 <- function(x,y,m,V,b){
    fxn <- y ~Co*((1/(1 + (m/V)*b)) + (((m/V)*b)/(1 + (m/V)*b)) * exp(-((1 + (m/V)*b)/((m/V)*b))*kfasS*x))
    grd1 <- data.frame(Co = c(0,1000),
                       kfasS  = c(0,10))
    cc<- capture.output(type="message",
                        fit9 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit9)==TRUE){
      fit9 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsfas <- as.vector(coefficients(fit9))
      pars_Co <- parsfas[1L]; pars_kfasS <- parsfas[2L]; pars_lin <- parsfas[3L]
      fun.1 <- function(x) {pars_Co*pars_lin*((1/(1 + (m/V)*b)) + (((m/V)*b)/(1 + (m/V)*b)) * exp(-((1 + (m/V)*b)/((m/V)*b))*pars_kfasS*x))}
      r <- fun.1(0)
      Comin = r*0.9                ;Comax = r*1.1
      kfasSmin =  pars_kfasS*0.9   ;kfasSmax =  pars_kfasS*1.1
      grd2 <- data.frame(Co = c(Comin,Comax),
                         kfasS = c(kfasSmin,kfasSmax))
      fit9 <- nls2(fxn,
                   start = grd2,
                   algorithm = "brute-force",
                   control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit9)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <- rmse(y,predict(fit9))
      mae  <- mae(y,predict(fit9))
      mse  <- mse(y,predict(fit9))
      rae  <- rae(y,predict(fit9))
      PAIC <- AIC(fit9)
      PBIC <- BIC(fit9)
      SE   <- sqrt((sum((y-predict(fit9))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Furusawa and Smith Model")
    print(summary(fit9))
    predval(x,n.dat)
    error(y)
    parsfas <- as.vector(coefficients(fit9))
    pars_Co <- parsfas[1L]; pars_kfasS <- parsfas[2L]
    theme_set(theme_bw())  # pre-set the bw theme.
    fun.1 <- function(x){pars_Co*((1/(1 + (m/V)*b)) + (((m/V)*b)/(1 + (m/V)*b)) * exp(-((1 + (m/V)*b)/((m/V)*b))*pars_kfasS*x))}
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of Ct vs time with non-linear Furusawa and Smith model",
           y="Ct",
           x="time",
           title="Furusawa and Smith Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
   }
  EQ2 <- function(x,y,m,V,b,Co){
    fxn <- y ~Co*((1/(1 + (m/V)*b)) + (((m/V)*b)/(1 + (m/V)*b)) * exp(-((1 + (m/V)*b)/((m/V)*b))*kfasS*x))
    grd1 <- data.frame(kfasS  = c(0,10))
    cc<- capture.output(type="message",
                        fit9 <- try(nls2::nls2(fxn,
                                               data = dat,
                                               start = grd1,
                                               algorithm = "port",
                                               control = list(maxiter = 1000)),
                                    silent=TRUE))
    if(is.null(fit9)==TRUE){
      fit9 <- nls2(fxn,
                   data = dat,
                   start = grd1,
                   algorithm = "plinear-random",
                   control = list(maxiter = 1000))
      parsfas <- as.vector(coefficients(fit9))
      pars_kfasS <- parsfas[1L]; pars_lin <- parsfas[2L]
      kfasSmin <-pars_kfasS*0.9   ; kfasSmax <-pars_kfasS*1.1
      grd2 <- data.frame(kfasS = c(kfasSmin,kfasSmax))
      fit9 <- nls2(fxn,
                     start = grd2,
                     algorithm = "brute-force",
                     control=list(maxiter=1000))
    }else{}
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.val <- predict(fit9)
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
      message("Estimated Values")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse <-rmse(y,predict(fit9))
      mae  <- mae(y,predict(fit9))
      mse  <- mse(y,predict(fit9))
      rae  <- rae(y,predict(fit9))
      PAIC <- AIC(fit9)
      PBIC <- BIC(fit9)
      SE   <- sqrt((sum((y-predict(fit9))^2))/(n.dat-2))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    message("Furusawa and Smith Model")
    print(summary(fit9))
    predval(x,n.dat)
    error(y)
    parsfas <- as.vector(coefficients(fit9))
    pars_kfasS <- parsfas[1L]
    theme_set(theme_bw())
    fun.1 <- function(x){Co*((1/(1 + (m/V)*b)) + (((m/V)*b)/(1 + (m/V)*b)) * exp(-((1 + (m/V)*b)/((m/V)*b))*pars_kfasS*x))}
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_function(color="red", fun=fun.1,size=1)+
      geom_point()+
      labs(subtitle="Plot of Ct vs time with non-linear Furusawa and Smith model",
           y="Ct",
           x="time",
           title="Furusawa and Smith Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  if(missing(Co)){
    EQ1(x,y,m,V,b)}
  else if(is.null(Co)){
    EQ1(x,y,m,V,b)}
  else if(isFALSE(Co)){
    EQ1(x,y,m,V,b)}
  else{EQ2(x,y,m,V,b,Co)
    }
}
