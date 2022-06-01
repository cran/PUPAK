#' @title Non-Linear Avrami Adsorption Kinetic Model
#' @description The Avrami Adsorption Kinetic Model investigates the time-concentration profiles of sorbent-sorbate interactions in adsorption-based water treatment. This equation was developed with the experimentally supported assumptions that the new phase is nucleated by germ nuclei that already exist in the old phase, and whose number can be altered by previous treatment (Oladoja, 2016).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qe the numerical value for the amount adsorbed at equilibrium
#' @import nls2
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the non-linear regression and the parameter estimation for the Avrami adsorption kinetic model analysis
#' @examples
#' \donttest{
#' t <- c(0,15,30,45,60,75,90,105,120)
#' qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qe <- 4.68
#' avrami(t,qt,qe)}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Lopes, E. C. N., dos Anjos, F. S. C., Vieira, E. F. S., &; Cestari, A. R. (2003) <doi:10.1016/S0021-9797(03)00326-6> An alternative Avrami equation to evaluate kinetic parameters of the interaction of Hg(II) with thin chitosan membranes. Journal of Colloid and Interface Science, 263(2), 542-547.
#' @references Oladoja, N. A. (2016) <doi:10.1080/19443994.2015.1076355> A critical review of the applicability of Avrami fractional kinetic equation in adsorption-based water treatment studies. Desalination and Water Treatment, 57(34), 15813-15825.
#' @export

avrami <- function(t,qt,qe){
  x     <- t
  y     <- qt
  dat   <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  fxn  <- y ~ (qe* (1- ((exp(-k1*x))^n)))
  grd1 <- data.frame(k1 = c(0,10),
                     n = c(0,1))
  cc<- capture.output(type="message",
                      fit4 <- try(nls2::nls2(fxn,
                                             data = dat,
                                             start = grd1,
                                             algorithm = "port",
                                             control = list(maxiter = 1000)),
                                  silent=TRUE))
  if(is.null(fit4)==TRUE){
    fit4 <- nls2(fxn,
                 data = dat,
                 start = grd1,
                 algorithm = "plinear-random",
                 control = list(maxiter = 1500))
    parsavm <- as.vector(coefficients(fit4))
    pars_k1 <- parsavm[1L]
    k1min <- pars_k1*0.9; k1max <- pars_k1*1.1
    grd2 <- data.frame(k1=c(k1min,k1max),
                       n = c(0,1))
    fit4 <- nls2(fxn,
                 start = grd2,
                 algorithm = "brute-force",
                 control=list(maxiter=1000))
  }else{}
  predval <- function(x,n.dat){
    Col1  <- c(rep(" |",each = n.dat))
    Col2  <- c(rep("|",each = n.dat))
    pred.val <- predict(fit4)
    time  <- x
    P.Table <- data.frame(Col1,time,Col1,pred.val,Col2)
    colnames(P.Table) <- c(" |","Time "," |","Pred Val","|")
    message("Estimated Values")
    print(P.Table, right=T, row.names = F)
  }
  error <- function(y){
    rmse <- rmse(y,predict(fit4))
    mae  <- mae(y,predict(fit4))
    mse  <- mse(y,predict(fit4))
    rae  <- rae(y,predict(fit4))
    PAIC <- AIC(fit4)
    PBIC <- BIC(fit4)
    SE   <- sqrt((sum((y-predict(fit4))^2))/(n.dat-2))
    Col1 <- c(" |"," |"," |"," |"," |"," |"," |")
    Col2 <- c("|","|","|","|","|","|","|")
    E.P  <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ")
    E.V  <- c(rmse,mae,mse,rae,PAIC,PBIC,SE)
    E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
    colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
    message("Error Estimation")
    print(E.Table, right=F, row.names = F)
  }
  message("Avrami Model")
  print(summary(fit4))
  predval(x,n.dat)
  error(y)
  parsavm1 <- as.vector(coefficients(fit4))
  pars_k1  <- parsavm1[1L]; pars_n <- parsavm1[2L]
  theme_set(theme_bw())
  fun.1 <- function(x) (qe* (1- (exp(-pars_k1*x)^pars_n)))
  plot <- ggplot(dat, aes(x=x,y=y))+
    geom_function(color="red", fun=fun.1, size=1)+
    geom_point()+
    labs(subtitle="Plot of qt vs time with non-linear Avrami Model",
         y="qt",
         x="time",
         title="Avrami Model",
         caption="Created by PUPAK using ggplot2")
  print(plot)
}
