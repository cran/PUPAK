#' @title Linear Boyd External Diffusion Model
#' @description  The Boyd External Diffusion Model is frequently applied to adsorption kinetic data to calculate the rate constant, assuming that film diffusion is the rate-limiting step in the first few minutes of the adsorption process. The film diffusion has a strong dependency on agitation. Boyd’s diffusion models are used in numerous adsorption studies mostly to determine the rate-controlling step (Viegas, Campinas, Costa, and Rosa, 2014).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t. 
#' @param qinf the numerical value for the amount adsorbed at infinite time. If the value for qinf is not given, the highest qt value will be considered as the qinf. 
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the linear regression and the parameter estimation for the Boyd External Diffusion Model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <-c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples qinf <- 4.68
#' @examples bed.lm(t,qt,qinf)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Boyd, G. E., Adamson, A. W., & Myers, L. S. (1947) <doi:10.1021/ja01203a066> The Exchange Adsorption of Ions from Aqueous Solutions by Organic Zeolites. II. Kinetics1. Journal of the American Chemical Society, 69(11), 2836-2848.
#' @references Viegas, R. M. C., Campinas, M., Costa, H., &; Rosa, M. J. (2014) <doi:10.1007/s10450-014-9617-9> How do the HSDM and Boyd's model compare for estimating intraparticle diffusion coefficients in adsorption processes. Adsorption, 20(5-6), 737-746.
#' @export

bed.lm<- function(t,qt,qinf){
  if (missing(qinf)){
    qinf <- max(qt)
  } 
  else if(is.null(qinf)){
    qinf <- max(qt)
  }
  else{qinf <- qinf}
  dat1 <- data.frame(t,qt)
  qt <- dat1$qt[which(dat1$qt < qinf)]
  x  <- dat1$t[which(dat1$qt < qinf)]
  y     <- log(1-(qt/qinf))
  dat1  <- data.frame(x,y)
  dat   <- subset(dat1,  y!="-Inf")
  dat   <- subset(dat,  y!="NaN")
  n1    <- nrow(dat1)
  n.dat <- nrow(dat)
  y     <- dat$y        ;x <- dat$x

  fit20   <- lm(y ~ x)
  lin.val <- coef(fit20)
  int <-as.numeric(lin.val[1])
  slp <-as.numeric(lin.val[2])

  R <- -slp
  A <- int

  pred.val <- (int+(slp*x))
  predval <- function(x,n.dat){
    Col1<- c(rep(" |",each = n.dat))
    Col2<- c(rep("|",each = n.dat))
    pred.qt <- qinf*((-exp(pred.val))+1)
    time <- x
    P.Table <- data.frame(Col1,x,Col1,pred.qt,Col1,pred.val,Col2)
    colnames(P.Table) <- c(" |","Time "," |","qt"," |","ln(1-(qt/qinf))","|")
    message("Estimated Values")
    print(P.Table, right=F, row.names = F)
  }

  error <- function(y){
    rmse   <- as.numeric(rmse(y,pred.val))
    mae    <- as.numeric(mae(y,pred.val))
    mse    <- as.numeric(mse(y,pred.val))
    rae    <- as.numeric(rae(y,pred.val))
    PAIC   <- as.numeric(AIC(fit20))
    PBIC   <- as.numeric(BIC(fit20))
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
  params <- function(R,A){
    param.name <- c("A=")
    param.val <- c(A)
    param.table <- data.frame(param.name,param.val)
    colnames(param.table) <- c("R=",R)
    message("BED Parameters")
    print(param.table, right=TRUE, row.names = F)
  }
  message("Boyd External Diffusion Model")
  message("Formula: qt = ln(1-(qt/qinf))=-Rt+A")
  message("Linear Model Summary")
  print(summary(fit20))
  params(R,A)
  predval(x,n.dat)
  error(y)

  xval <- seq(min(x),max(x),length=10000)
  eqm <- slp*xval+ int
  plot.dat <- data.frame (xval,eqm)

  theme_set(theme_bw())
  plot <- ggplot(dat, aes(x=x,y=y))+
    geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red") +
    geom_point(size=2)+
    labs(subtitle=expression("Plot of"~plain(ln)*bgroup("(",1-{frac(q[t],q[inf])},")")~"vs t with the linear Boyd External Diffusion model"),
         y=expression(plain(ln)*bgroup("(",1-{frac(q[t],q[inf])},")")),
         x="time",
         title="Boyd External Diffusion Linear Model",
         caption="Created by PUPAK using ggplot2")
  print(plot)
}
