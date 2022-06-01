#' @title Linear Elovich Adsorption Kinetic Model
#' @description The Elovich Adsorption Kinetic Model is an empirical rate equation that states the adsorption energy rises in a linear relationship with surface coverage. The model assumes that adsorption occurs on localized sites, the interaction between adsorbed ions is present, and the concentration of adsorbate is considered to be constant. It is applicable in gas adsorptions as well as wastewater processes (Largitte and Pasquier, 2016).
#' @param t the numerical value for contact time. This parameter should not be equal to zero to prevent infinite value. Any row(s) that contain(s) value of t equal to zero will be automatically removed to proceed with the calculation.
#' @param qt the numerical value for the amount adsorbed at time t
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the linear regression and the parameter estimation for the Elovich Adsorption Kinetic Model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <-c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples elovich.lm(t,qt)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Elovich, S. Y., & Larinov, O. G. (1962) Theory of adsorption from solutions of non electrolytes on solid (I) equation adsorption from solutions and the analysis of its simplest form,(II) verification of the equation of adsorption isotherm from solutions. Izvestiya Akademii Nauk, 2(2), 209–216.
#' @references Largitte, L., &; Pasquier, R. (2016) <doi:10.1016/j.cherd.2016.02.006> A review of the kinetics adsorption models and their application to the adsorption of lead by an activated carbon. Chemical Engineering Research and Design, 109, 495-504.
#' @export
elovich.lm<- function(t,qt){
  x <- log(t)   ;y <- qt
  dat1  <- data.frame(x,y)
  dat   <- subset(dat1, x!="-Inf")
  n1     <- nrow(dat1)
  n.dat <- nrow(dat)
  y <- dat$y    ;x <- dat$x

  fit18 <- lm(y ~ x)
  lin.val <- coef(fit18)
  int<-as.numeric(lin.val[1])
  slp <-as.numeric(lin.val[2])
  b <- 1/slp
  a <- exp(int*b) + b

  pred.val <- (int+(slp*x))
  predval <- function(x,n.dat){
    Col1<- c(rep(" |",each = n.dat))
    Col2<- c(rep("|",each = n.dat))
    pred.qt <- pred.val
    time <- exp(x)
    P.Table <- data.frame(Col1,time,Col1,x,Col1,pred.qt,Col2)
    colnames(P.Table) <- c(" |","Time "," |","ln(Time)"," |","qt","|")
    message("Estimated Values")
    print(P.Table, right=F, row.names = F)
  }
  error <- function(y){
    rmse   <- as.numeric(rmse(y,pred.val))
    mae    <- as.numeric(mae(y,pred.val))
    mse    <- as.numeric(mse(y,pred.val))
    rae    <- as.numeric(rae(y,pred.val))
    PAIC   <- as.numeric(AIC(fit18))
    PBIC   <- as.numeric(BIC(fit18))
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
    message("Elovich Parameters")
    print(param.table, right=TRUE, row.names = F)
  }
  message("Elovich Model")
  message("Formula: qt = (ln(ab)/b)+(ln(t)/b) Assuming abt >> 1")
  message("Linear Model Summary")
  print(summary(fit18))
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
    labs(subtitle="Plot of the qt vs ln(t) with linear Elovich model",
         y="qt",
         x="ln(time)",
         title="Elovich Linear Model",
         caption="Created by PUPAK using ggplot2")
  print(plot)
}
