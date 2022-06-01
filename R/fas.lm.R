#' @title Linear Furusawa and Smith Diffusion Model
#' @description The Furusawa and Smith Diffusion Model is known to describe the rate of adsorption, assuming that only external diffusion resistance was predominant during the initial adsorption period and controlled the adsorption rate. The diffusion model relates the change in fluid phase concentration and time with the fluid phase concentration at the external surface and an external mass transfer coefficient (Furusawa & Smith, 1974).
#' @param t the numerical value for contact time
#' @param Ct the numerical value for the concentration of the adsorbent at time t. This parameter should not contain a value equal to zero. Any row(s) that contain(s) value of Ct equal to zero will be automatically removed to proceed with the calculation.
#' @param m the numerical value for mass of adsorbent
#' @param V the numerical value for volume of solution
#' @param b the numerical value for the Langmuir isotherm constant
#' @param Co the numerical value for the initial concentration of the adsorbent
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @return the linear regression and the parameter estimation for the Furusawa and Smith Diffusion Model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples Ct <- c(10.000,8.141,8.056,7.949,7.863,7.799,7.778,7.756,7.692)
#' @examples m <- 0.05
#' @examples V <- 0.1
#' @examples b <- 1.3
#' @examples Co <- 10
#' @examples fas.lm(t,Ct,m,V,b,Co)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Furusawa, T., & Smith, J. M. (1974) <doi:10.1002/aic.690200111> Intraparticle mass transport in slurries by dynamic adsorption studies. AIChE Journal, 20(1), 88–93.
#' @references Yakub, E., Agarry, S. E., Omoruwou, F., &; Owabor, C. N. (2020) <doi:10.1080/02726351.2019.1616862> Comparative study of the batch adsorption kinetics and mass transfer in phenol-sand and phenol-clay adsorption systems. Particulate Science and Technology, 38(7).
#' @export
fas.lm<- function(t,Ct,m,V,b,Co){
  m <- m   ;V  <- V   ;b <- b
  x <- t   ;Co <- Co
  y <- (1/(1+(1/((m/V)*b))))*log((Ct/Co)-((1/((m/V)*b))*(1-(Ct/Co))))
  dat1 <- data.frame(x,y)
  dat   <- subset(dat1,  y!="-Inf")
  dat   <- subset(dat,  y!="NaN")
  n1    <- nrow(dat1)
  n.dat <- nrow(dat)
  y <- dat$y  ;x <- dat$x

  fit21  <- lm(y ~ x)
  lin.val <- coef(fit21)
  int <- as.numeric(lin.val[1])
  slp <- as.numeric(lin.val[2])
  kFaSS <- -slp
  pred.val <- (int+(slp*x))

  predval  <- function(x,n.dat){
    Col1 <- c(rep(" |",each = n.dat))
    Col2 <- c(rep("|",each  = n.dat))
    pred.Ct <- Co*((1/(1 + (m/V)*b)) + (((m/V)*b)/(1 + (m/V)*b)) * exp(-((1 + (m/V)*b)/((m/V)*b))*kFaSS*x))
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
    PAIC   <- as.numeric(AIC(fit21))
    PBIC   <- as.numeric(BIC(fit21))
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
  params <- function(kFaSS){
    param.name <- c("")
    param.val <- c("")
    param.table <- data.frame(param.name,param.val)
    colnames(param.table) <- c("kFaSS=",kFaSS)
    message("FaS Parameters")
    print(param.table, right=TRUE, row.names = F)
  }
  message("Furusawa and Smith Diffusion Model")
  message("Formula: (1/(1+(1/((m/V)*b))))*log((Ct/Co)-((1/((m/V)*b))*(1-(Ct/Co))))=-kFaSS*t")
  message("Linear Model Summary")
  print(summary(fit21))
  params(kFaSS)
  predval(x,n.dat)
  error(y)

  xval <- seq(min(x),max(x),length=10000)
  eqm  <- slp*xval+ int
  plot.dat <- data.frame (xval,eqm)

  theme_set(theme_bw())
  plot <- ggplot(dat, aes(x=x,y=y))+
    geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red") +
    geom_point(size=2)+
    labs(subtitle=expression("Plot of"~frac(1,1+frac(1,frac(m,v)*b))*plain(ln)*
                             bgroup("(",{frac(C[t],C[o])-frac(1,frac(m,v)*b)*
                                 bgroup("(",1-frac(C[t],C[o]),")")},")")~"vs t with the linear Furasawa and Smith Diffusion model"),
         y=expression(paste(frac(1,1+frac(1,frac(m,v)*b))*plain(ln)*
                            bgroup("(",{frac(C[t],C[o])-frac(1,frac(m,v)*b)*
                                bgroup("(",1-frac(C[t],C[o]),")")},")"))),
         x="time",
         title="Furasawa and Smith Diffusion Linear Model",
         caption="Created by PUPAK using ggplot2")
  print(plot)
}
