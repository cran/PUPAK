#' @title Linear Matthews and Weber Diffusion Model
#' @description The Matthews and Weber Diffusion Model is used in determining the rate constant for film diffusion where it assumes that intraparticle diffusion can be neglected at the early period of contact. Derived from Fickien’s law application, the solute concentration in the liquid phase in this model is expressed as a function of solute concentration difference in the liquid phase and at the adsorbent surface (Prasad & Srivastava, 2009).
#' @param t the numerical value for contact time
#' @param Ct the numerical value for the concentration of the adsorbent at time t. This parameter should not contain a value equal to zero. Any row(s) that contain(s) value of Ct equal to zero will be automatically removed to proceed with the calculation.
#' @param Co the numerical value for the initial concentration of the adsorbent
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the linear regression and the parameter estimation for the Matthews and Weber diffusion model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples Ct <-c(10.000,8.141,8.056,7.949,7.863,7.799,7.778,7.756,7.692)
#' @examples Co <- 10
#' @examples maw.lm(t,Ct,Co)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Mathews, A. P., &; Weber, W. J. (1984) <doi:10.1080/00986448408940104> Modeling and parameter evaluation for adsorption in slurry reactors. Chemical Engineering Communications, 25(1-6), 157-171.
#' @references Krishna Prasad, R., & Srivastava, S. N. (2009) <doi:10.1016/j.cej.2008.05.021> Sorption of distillery spent wash onto fly ash: Kinetics and mass transfer studies. Chemical Engineering Journal, 146(1), 90–97.
#' @export

maw.lm<- function(t,Ct,Co){
  Co <- Co        ;x <- t
  y  <- log(Ct/Co)
  dat1  <- data.frame(x,y)
  dat   <- subset(dat1,  y!="-Inf")
  dat   <- subset(dat,  y!="NaN")
  n1    <- nrow(dat1)
  n.dat <- nrow(dat)
  y <- dat$y      ;x <- dat$x

  fit22  <- lm(y ~ x)
  lin.val <- coef(fit22)
  int <-as.numeric(lin.val[1])
  slp <-as.numeric(lin.val[2])
  kMaWS   <- -slp
  pred.val <- (int+(slp*x))

  predval <- function(x,n.dat){
    Col1 <- c(rep(" |",each = n.dat))
    Col2 <- c(rep("|",each = n.dat))
    pred.Ct <- exp(pred.val)*Co
    time <- x
    P.Table <- data.frame(Col1,x,Col1,pred.Ct,Col1,pred.val,Col2)
    colnames(P.Table) <- c(" |","Time "," |","Ct"," |","ln(Ct/Co)","|")
    message("Estimated Values")
    print(P.Table, right=F, row.names = F)
  }

  error <- function(y){
    rmse   <- as.numeric(rmse(y,pred.val))
    mae    <- as.numeric(mae(y,pred.val))
    mse    <- as.numeric(mse(y,pred.val))
    rae    <- as.numeric(rae(y,pred.val))
    PAIC   <- as.numeric(AIC(fit22))
    PBIC   <- as.numeric(BIC(fit22))
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
  params <- function(kMaWS){
    param.name <- c("")
    param.val <- c("")
    param.table <- data.frame(param.name,param.val)
    colnames(param.table) <- c("kMaWS=",kMaWS)
    message("MaW Parameter")
    print(param.table, right=TRUE, row.names = F)
    }

  message("Matthews and Weber Diffusion Model")
  message("Formula: ln(Ct/Co)=-kMaWS*t")
  message("Linear Model Summary")
  print(summary(fit22))
  params(kMaWS)
  predval(x,n.dat)
  error(y)

  xval <- seq(min(x),max(x),length=10000)
  eqm = slp*xval+ int
  plot.dat <- data.frame (xval,eqm)

  theme_set(theme_bw())
  plot <- ggplot(dat, aes(x=x,y=y))+
    geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red") +
    geom_point(size=2)+
    labs(subtitle=expression("Plot of"~plain(ln)*bgroup("(",1-{frac(C[t],C[0])},")")~"vs t with the linear Matthews and Weber Diffusion model"),
         y=expression(plain(ln)*bgroup("(",1-{frac(C[t],C[0])},")")),
         x="time",
         title="Matthews and Weber Diffusion Linear Model",
         caption="Created by PUPAK using ggplot2")
  print(plot)
}
