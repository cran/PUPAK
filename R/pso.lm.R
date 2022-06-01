#' @title Linear Pseudo-Second-Order Adsorption Kinetic Model
#' @description The Pseudo-Second-Order Adsorption Kinetic Model is an empirical rate equation known to be the simplified second-order expression of the Pseudo-First Order Adsorption Kinetic Model. It is widely applied to adsorption systems, from biomass to nanomaterials as adsorbent and from heavy metals to pharmaceuticals as adsorbate or contaminant (Revellame, Fortela, Sharp, Hernandez, and Zappi, 2020).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t. This parameter should not contain a value equal to zero as it will cause an infinite value. Any row(s) that contain(s) value of qt equal to zero will be automatically removed to proceed with the calculation.
#' @param qe the numerical value for the amount adsorbed at equilibrium
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the linear regression and the parameter estimation for the Pseudo-Second-Order Model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <-c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples qe <- 4.8
#' @examples pso.lm(t,qt,qe)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Ho, Y. S., &; Mckay, G. (1999) <doi:10.1016/S0032-9592(98)00112-5> Pseudo-second order model for sorption processes. In Process Biochemistry (Vol. 34).
#' @references Revellame, E. D., Fortela, D. L., Sharp, W., Hernandez, R., &; Zappi, M. E. (2020) <doi:10.1016/j.clet.2020.100032>. Adsorption kinetic modeling using pseudo-first order and pseudo-second order rate laws: A review. In Cleaner Engineering and Technology (Vol. 1). Elsevier Ltd.
#' @export
pso.lm<- function(t,qt,qe){
  x <- t          ;y <- t/qt
  dat1  <- data.frame(x,y)
  dat   <- dat1[complete.cases(dat1),]
  n1    <- nrow(dat1)
  n.dat <- nrow(dat)
  y <- dat$y      ;x <- dat$x
  EQ1 <- function(x,y){
    fit17   <- lm(y ~ x)
    lin.val <- coef(fit17)
    int <-as.numeric(lin.val[1])
    slp <-as.numeric(lin.val[2])
    qe  <-(1/slp)
    k2  <-(1/((qe^2)*int))
    pred.val <- ((slp*x)+int)

    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      pred.qt <- x/pred.val
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.qt,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","qt"," |","time/qt","|")
      message("Estimated Values")
      print(P.Table, right=F, row.names = F)
    }

    error <- function(y){
      rmse   <- as.numeric(rmse(y,pred.val))
      mae    <- as.numeric(mae(y,pred.val))
      mse    <- as.numeric(mse(y,pred.val))
      rae    <- as.numeric(rae(y,pred.val))
      PAIC   <- as.numeric(AIC(fit17))
      PBIC   <- as.numeric(BIC(fit17))
      SE     <- as.numeric(sqrt((sum((y-pred.val)^2))/(n.dat-2)))
      rsqtot <- as.numeric(cor(x,y)^2)
      Col1  <- c(" |"," |"," |"," |"," |"," |"," |"," |")
      Col2  <- c("|","|","|","|","|","|","|","|")
      E.P   <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ","Coefficient of Determination (R^2) ")
      E.V   <- c(rmse,mae,mse,rae,PAIC,PBIC,SE,rsqtot)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    params <- function(qe,k2){
      param.name <- c("k2=")
      param.val <- c(k2)
      param.table <- data.frame(param.name,param.val)
      colnames(param.table) <- c("qe=",qe)
      message("PSO Parameters")
      print(param.table, right=TRUE, row.names = F)
    }

    message("Pseudo-Second-Order Model")
    message("Formula: (t/qt)=(1/(qe^2*k2))+(t/qe)")
    message("Linear Model Summary")
    print(summary(fit17))
    params(qe,k2)
    predval(x,n.dat)
    error(y)

    xval <- seq(min(x),max(x),length=10000)
    eqm <- slp*xval+ int
    plot.dat <- data.frame (xval,eqm)

    theme_set(theme_bw())
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red")+
      geom_point(size=2)+
      labs(subtitle="Plot of the t/qt vs t with linear Pseudo Second Order model",
           y="time/qt",
           x="time",
           title="Pseudo Second Order Linear Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
  EQ2 <- function(x,y,qe){
    slp <- (1/qe)
    fit17 <- lm(I(y - slp*x) ~ 1)
    lin.val  <- coef(fit17)
    int <-as.numeric(lin.val[1])
    pred.val <- ((slp*x)+int)
    RSS <- sum((y-pred.val)^2)
    TSS <- sum((y-mean(y))^2)

    qe <- (1/slp)
    k2 <- (1/((qe^2)*int))

    predval <- function(x,n.dat){
      Col1<- c(rep(" |",each = n.dat))
      Col2<- c(rep("|",each = n.dat))
      pred.qt <- x/pred.val
      time <- x
      P.Table <- data.frame(Col1,time,Col1,pred.qt,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","qt"," |","time/qt","|")
      message("Estimated Values")
      print(P.Table, right=F, row.names = F)
    }

    error <- function(y){
      rmse   <- as.numeric(rmse(y,pred.val))
      mae    <- as.numeric(mae(y,pred.val))
      mse    <- as.numeric(mse(y,pred.val))
      rae    <- as.numeric(rae(y,pred.val))
      PAIC   <- as.numeric(AIC(fit17))
      PBIC   <- as.numeric(BIC(fit17))
      SE     <- as.numeric(sqrt((sum((y-pred.val)^2))/(n.dat-2)))
      rsqtot <- as.numeric(1-(RSS/TSS))
      Col1 <- c(" |"," |"," |"," |"," |"," |"," |"," |")
      Col2 <- c("|","|","|","|","|","|","|","|")
      E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ","Coefficient of Determination (R^2) ")
      E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE,rsqtot)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation")
      print(E.Table, right=F, row.names = F)
    }
    params <- function(qe,k2){
      param.name <- c("k2=")
      param.val <- c(k2)
      param.table <- data.frame(param.name,param.val)
      colnames(param.table) <- c("qe=",qe)
      message("PSO Parameters")
      print(param.table, right=TRUE, row.names = F)
    }

    message("Pseudo-Second-Order Model")
    message("Formula: (t/qt)=(1/(qe^2*k2))+(t/qe)")
    message("Linear Model Summary")
    print(summary(fit17))
    params(qe,k2)
    predval(x,n.dat)
    error(y)

    xval <- seq(min(x),max(x),length=10000)
    eqm = slp*xval+ int
    plot.dat <- data.frame (xval,eqm)

    theme_set(theme_bw())
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red")+
      geom_point(size=2)+
      labs(subtitle="Plot of the t/qt vs t with linear Pseudo-Second-Order model",
           y="time/qt",
           x="time",
           title="Pseudo Second Order Linear Model",
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
