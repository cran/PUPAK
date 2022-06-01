#' @title Linear Pseudo-First-Order Adsorption Kinetic Model
#' @description The Pseudo-First-Order Adsorption Kinetic Model follows the Linear Driving Force model (LDF) which states that the rate of mass transfer is equal to the transfer coefficient and the difference between the amount adsorbed and the amount adsorbed in equilibrium. This model is an empirical rate equation known to describe the rate of sorption in liquid-phase systems. The PFO model is not suitable for the whole adsorption reaction since the rate of adsorption decreases until it reaches the maximum adsorption capacity, and thus, the rate is zero at equilibrium (Plazinski, Rudzinski, and Plazinska, 2009).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t. This parameter should not be greater than or equal to qe as it will cause an incalculable value. Any row(s) that contain(s) a value greater than or equal to qe will be automatically removed to proceed with the calculation.
#' @param qe the numerical value for the amount adsorbed at equilibrium
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the linear regression and the parameter estimation for the Pseudo First Order Model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <-c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples qe <- 4.8
#' @examples pfo.lm(t,qt,qe)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Lagergren, S. (1898), Zur theorie der sogenannten adsorption gelster stoffe, Kungliga Svenska Vetenskapsakademiens. Handlingar, 24 (4) : 1-39.
#' @references Plazinski, W., Rudzinski, W., &; Plazinska, A. (2009) <doi:10.1016/j.cis.2009.07.009> Theoretical models of sorption kinetics including a surface reaction mechanism: A review. In Advances in Colloid and Interface Science</i> (Vol. 152, Issues 1-2, pp. 2-13).
#' @export
pfo.lm<- function(t,qt,qe){
  qe <- qe         ;x   <- t
  y  <- log(qe-qt) ;int <- log(qe)
  dat1  <- data.frame(x,y)
  dat   <- dat1[complete.cases(dat1),]
  n1    <- nrow(dat1)
  n.dat <- nrow(dat)
  y <- dat$y       ;x <- dat$x
  fit16   <- lm(I(y -int) ~ 0+x)
  lin.val <- coef(fit16)
  slp <-as.numeric(lin.val[1])
  k1  <- -slp
  pred.val <- (int+(slp*x))
  predval <- function(x,n.dat){
    Col1 <- c(rep(" |",each = n.dat))
    Col2 <- c(rep("|",each = n.dat))
    pred.qt <- -exp(pred.val)+qe
    time <- x
    P.Table <- data.frame(Col1,time,Col1,pred.qt,Col1,pred.val,Col2)
    colnames(P.Table) <- c(" |","Time "," |","qt"," |","ln(qe-qt)","|")
    message("Estimated Values")
    print(P.Table, right=F, row.names = F)
  }
  error <- function(y){
    rmse   <- as.numeric(rmse(y,pred.val))
    mae    <- as.numeric(mae(y,pred.val))
    mse    <- as.numeric(mse(y,pred.val))
    rae    <- as.numeric(rae(y,pred.val))
    PAIC   <- as.numeric(AIC(fit16))
    PBIC   <- as.numeric(BIC(fit16))
    SE     <- as.numeric(sqrt((sum((y-pred.val)^2))/(n.dat-2)))
    rsqtot <- as.numeric(summary(fit16)$r.squared)
    Col1 <- c(" |"," |"," |"," |"," |"," |"," |"," |")
    Col2 <- c("|","|","|","|","|","|","|","|")
    E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ","Coefficient of Determination (R^2) ")
    E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE,rsqtot)
    E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
    colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
    message("Error Estimation")
    print(E.Table, right=F, row.names = F)
  }
  params <- function(k1,qe){
    param.name <- c("k1=")
    param.val <- c(k1)
    param.table <- data.frame(param.name,param.val)
    colnames(param.table) <- c("qe=",qe)
    message("PFO Parameters")
    print(param.table, right=TRUE, row.names = F)
  }
  message("Pseudo-First-Order Model")
  message("Formula: ln(qe-qt)=ln(qe)-k1*t")
  message("Linear Model Summary")
  print(summary(fit16))
  params(k1,qe)
  predval(x,n.dat)
  error(y)
  xval <- seq(min(x),max(x),length=10000)
  eqm = slp*xval+ int
  plot.dat <- data.frame (xval,eqm)
  theme_set(theme_bw())
  plot <- ggplot(dat, aes(x=x,y=y))+
    geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red")+
    geom_point(size=2)+
    labs(subtitle="Plot of the ln(qe-qt) vs t with linear Pseudo-First-Order model",
         y="ln(qe-qt)",
         x="time",
         title="Pseudo-First-Order Linear Model",
         caption="Created by PUPAK using ggplot2")
  print(plot)
}
