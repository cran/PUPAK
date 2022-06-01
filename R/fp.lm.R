#' @title Linear Fractional Power Adsorption Kinetic Model
#' @description The Fractional Power Adsorption Kinetic Model is an empirical rate equation in which the specific adsorption rate at a unit time can be estimated using its product's constant. The model was first used to explain the adsorption kinetics of phosphate in soils by anion exchange resin. It is also called the modified Freundlich equation (Netzahuatl-Muñoz, Del Carmen Cristiani-Urbina, and Cristiani-Urbina, 2015).
#' @param t the numerical value for contact time. This parameter should not contain a value equal to zero to prevent infinite value. Any row(s) that contain(s) value of t equal to zero will be automatically removed to proceed with the calculation.
#' @param qt the numerical value for the amount adsorbed at time t. This parameter should not contain a value equal to zero to prevent incalculable value. Any row(s) that contain(s) value of qt equal to zero will be automatically removed to proceed with the calculation.
#' @param qe the numerical value for the amount adsorbed at equilibrium
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the linear regression and the parameter estimation for the Fractional Power model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <-c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples qe <- 4.8
#' @examples fp.lm(t,qt,qe)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Dalai, R. C. (1974) <doi: 10.1080/00103627409366531> Desorption of soil phosphate by anion-exchange resin. Communications in Soil Science and Plant Analysis, 5(6), 531-538.
#' @references Netzahuatl-Munoz, A. R., del Carmen Cristiani-Urbina, M., &; Cristiani-Urbina, E. (2015) <doi:10.1371/journal.pone.0137086> Chromium Biosorption from Cr(VI) aqueous solutions by Cupressus lusitanica bark: Kinetics, equilibrium and thermodynamic studies. PLoS ONE, 10(9).
#' @export
fp.lm<- function(t,qt,qe){
  qe <- qe       ;x <- log(t)
  y  <- log(qt/qe)
  dat1  <- data.frame(x,y)
  dat   <- subset(dat1, x!="-Inf")
  dat   <- subset(dat,  y!="-Inf")
  n1    <- nrow(na.omit(dat1))
  n.dat <- nrow(na.omit(dat))
  y  <- dat$y    ;x <- dat$x

  fit19   <- lm(y ~ x)
  lin.val <- coef(fit19)
  int <-as.numeric(lin.val[1])
  slp <-as.numeric(lin.val[2])
  b <- slp  ;a <- exp(int)

  pred.val <- (int+(slp*x))
  predval <- function(x,n.dat){
    Col1 <- c(rep(" |",each = n.dat))
    Col2 <- c(rep("|",each = n.dat))
    pred.qt <- exp(pred.val)*qe
    time <- exp(x)
    P.Table <- data.frame(Col1,time,Col1,x,Col1,pred.qt,Col1,pred.val,Col2)
    colnames(P.Table) <- c(" |","Time "," |","ln(Time)"," |","qt"," |","ln(qt/qe)","|")
    message("Estimated Values")
    print(P.Table, right=F, row.names = F)
  }
  error <- function(y){
    rmse   <- as.numeric(rmse(y,pred.val))
    mae    <- as.numeric(mae(y,pred.val))
    mse    <- as.numeric(mse(y,pred.val))
    rae    <- as.numeric(rae(y,pred.val))
    PAIC   <- as.numeric(AIC(fit19))
    PBIC   <- as.numeric(BIC(fit19))
    SE     <- as.numeric(sqrt((sum((y-pred.val)^2))/(n.dat-2)))
    rsqtot <- as.numeric(cor(x,y)^2)
    Col1 <- c(" |"," |"," |"," |"," |"," |"," |"," |")
    Col2 <- c("|","|","|","|","|","|","|","|")
    E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ","Coefficient of Determination (R^2) ")
    E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE,rsqtot)
    E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
    colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
    message("Error Estimation")
    print(E.Table, right=F, row.names = F)
  }
  params <- function(a,b){
    param.name <- c("b=")
    param.val <- c(b)
    param.table <- data.frame(param.name,param.val)
    colnames(param.table) <- c("a=",a)
    message("Fractional Power Model Parameters")
    print(param.table, right=TRUE, row.names = F)
    }

  message("Fractional Power Model")
  message("Formula: ln(qt/qe)=ln(a)+bln(t)")
  message("Linear Model Summary")
  print(summary(fit19))
  params(a,b)
  predval(x,n.dat)
  error(y)

  xval <- seq(min(x),max(x),length=10000)
  eqm = slp*xval+ int
  plot.dat <- data.frame (xval,eqm)

  theme_set(theme_bw())
  plot <- ggplot(dat, aes(x=x,y=y))+
    geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red") +
    geom_point(size=2)+
    labs(subtitle="Plot of the ln(qt/qe) vs ln(t) with linear Fractional Power model",
         y="ln(qt/qe)",
         x="ln(time)",
         title="Fractional Power Linear Model",
         caption="Created by PUPAK using ggplot2")
  print(plot)
}
