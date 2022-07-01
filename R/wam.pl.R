#' @title Piecewise-Linear Weber and Morris Intraparticle Diffusion Model
#' @description The Weber and Morris Intraparticle Diffusion Model is a multi-linear sorption process in which the intraparticle diffusion process is the limiting factor, where other interaction mechanisms such as adsorption on the external surface and diffusion into the interior may be happening simultaneously (Campos, Barbosa, Rodriguez-Diaz, and Duarte, 2018).
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @import segmented
#' @import stats
#' @import ggplot2
#' @import Metrics
#' @import utils
#' @return the piecewise-linear regression and the parameter estimation for the Weber and Morris Intraparticle Diffusion Model
#' @examples t <- c(0,15,30,45,60,75,90,105,120)
#' @examples qt <-c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' @examples wam.pl(t,qt)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. Dela Cruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @references Weber, W.J. and Morris, J.C. (1963) Kinetics of Adsorption on Carbon from Solutions. Journal of the Sanitary Engineering Division, American Society of Civil Engineers, 89, 31-60.
#' @references Campos, N. F., Barbosa, C. M. B. M., Rodriguez-Diaz, J. M., &; Duarte, M. M. M. B. (2018) <doi:10.1177/0263617418773844> Removal of naphthenic acids using activated charcoal: Kinetic and equilibrium studies. Adsorption Science and Technology, 36(7-8), 1405-1421.
#' @export
wam.pl<- function(t,qt){
  x     <- sqrt(t) ;y <- qt
  dat   <- data.frame(x,y)
  n.dat <- nrow(na.omit(dat))
  fit15 <- lm(y ~ x)
  fit16 <- as.character(segmented(fit15, seg.Z=~x, psi=x[3], control=seg.control(display=FALSE,)), K=3)
  if(length(fit16)==12){
    fit15    <- lm(y ~ x)
    lin.val  <- coef(fit15)
    int      <- as.numeric(lin.val[1])
    slp      <- as.numeric(lin.val[2])
    pred.val <- ((slp*x)+int)
    params <- function(slp,int){
      param.name <- c("kWAM=","C=")
      param.val <- c(slp,int)
      param.table <- data.frame(param.name,param.val)
      colnames(param.table) <- c("","")
      message("Weber and Morris Parameters")
      print(param.table, right=TRUE, row.names = F)
    }
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      time <- x^2
      sqtime <- x
      P.Table <- data.frame(Col1,time,Col1,sqtime,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","sqrt(Time)"," |","qt","|")
      message("Estimated Values:")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse    <- round(as.numeric(rmse(y,pred.val)),digits=10)
      mae     <- round(as.numeric(mae(y,pred.val)),digits=10)
      mse     <- round(as.numeric(mse(y,pred.val)),digits=10)
      rae     <- round(as.numeric(rae(y,pred.val)),digits=10)
      PAIC    <- round(as.numeric(AIC(fit15)),digits=10)
      PBIC    <- round(as.numeric(BIC(fit15)),digits=10)
      SE      <- round(as.numeric(sqrt((sum((y-pred.val)^2))/(n.dat-2))),digits=10)
      rsqtot  <- round(as.numeric(cor(x,y)^2),digits=10)
      Col1    <- c(" |"," |"," |"," |"," |"," |"," |"," |")
      Col2    <- c("|","|","|","|","|","|","|","|")
      E.P     <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ","Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ","Standard Error Estimate ","Coefficient of Determination (R^2) ")
      E.V     <- c(rmse,mae,mse,rae,PAIC,PBIC,SE,rsqtot)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation:")
      print(E.Table, right=F, row.names = F)
    }
    message("Weber-Morris Diffusion Model")
    message("Formula: qt =  (kWAM*sqrt(t))+C")
    message("Linear Model Summary")
    print(summary(fit15))
    params(slp,int)
    predval(x,n.dat)
    error(y)
    xval <- seq(min(x),max(x),length=100)
    eqm = slp*xval+ int
    plot.dat <- data.frame (xval,eqm)
    theme_set(theme_bw())
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_line(data=plot.dat, aes(xval,eqm),size=1, color="red")+
      geom_point(size=2)+
      labs(subtitle="Plot of the qt vs square root of time with linear Weber-Morris model",
           y="qt",
           x="sqrt(time)",
           title="Weber-Morris Intraparticle Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }else{
    fit16 <- segmented(fit15, seg.Z=~x, psi=x[3], control=seg.control(display=FALSE,), K=3)
    brk.pts   <- fit16$psi[, 2]
    seg.range <- fit16$rangeZ[, 1]
    b <- intercept(fit16) ;y.int <- unlist(b)
    a <- slope(fit16)     ;slp <- unlist(a)
    int1 <- y.int[1]       ;slp1 <- slp[1]
    int2 <- y.int[2]       ;slp2 <- slp[2]
    brk <- brk.pts[1]
    firstlinearx  <- x[which(x < brk)]
    secondlinearx <- x[which(x > brk)]
    firstlineary  <- y[which(dat$x < brk)]
    secondlineary <- y[which(dat$x > brk)]
    low.pred  <- ((slp1*firstlinearx)+int1)
    high.pred <- ((slp2*secondlinearx)+int2)
    pred.val  <- c(low.pred,high.pred)
    RSS1 <- sum((low.pred-mean(firstlineary))^2)
    RSS2 <- sum((high.pred-mean(secondlineary))^2)
    TSS1 <- sum((firstlineary-mean(firstlineary))^2)
    TSS2 <- sum((secondlineary-mean(secondlineary))^2)
    params1 <- function(slp1,int1,slp2,int2){
      param.name <- c("kWAM1=","kWAM2=","C1=","C2=","Breakpoint")
      param.val <- c(slp1,slp2,int1,int2,brk.pts)
      param.table <- data.frame(param.name,param.val)
      colnames(param.table) <- c("","")
      message("Weber and Morris Parameters")
      print(param.table, right=TRUE, row.names = F)
    }
    predval <- function(x,n.dat){
      Col1 <- c(rep(" |",each = n.dat))
      Col2 <- c(rep("|",each = n.dat))
      time <- x^2
      sqtime <- x
      P.Table <- data.frame(Col1,time,Col1,sqtime,Col1,pred.val,Col2)
      colnames(P.Table) <- c(" |","Time "," |","sqrt(Time)"," |","qt","|")
      message("Estimated Values:")
      print(P.Table, right=T, row.names = F)
    }
    error <- function(y){
      rmse   <- as.numeric(rmse(y,pred.val))
      mae    <- as.numeric(mae(y,pred.val))
      mse    <- as.numeric(mse(y,pred.val))
      rae    <- as.numeric(rae(y,pred.val))
      PAIC   <- as.numeric(AIC(fit16))
      PBIC   <- as.numeric(BIC(fit16))
      SE     <- as.numeric(sqrt((sum((y-pred.val)^2))/(n.dat-2)))
      rsq1   <- as.numeric(cor(firstlinearx,firstlineary)^2)
      rsq2   <- as.numeric(cor(secondlinearx,secondlineary)^2)
      rsqtot <- as.numeric(summary(fit16)$r.squared)
      Col1   <- c(" |"," |"," |"," |"," |"," |"," |"," |"," |"," |")
      Col2   <- c("|","|","|","|","|","|","|","|","|","|")
      E.P <- c("Relative Mean Square Error ", "Mean Absolute Error ","Mean Squared Error ",
               "Relative Absolute Error ","Akaike Information Criterion ","Bayesian Information Criterion ",
               "Standard Error Estimate ","First Coefficient of Determination (R^2) ",
               "Second Coefficient of Determination (R^2) ","Overall Coefficient of Determination (R^2) ")
      E.V <- c(rmse,mae,mse,rae,PAIC,PBIC,SE,rsq1,rsq2,rsqtot)
      E.Table <- data.frame(Col1,E.P,Col1,E.V,Col2)
      colnames(E.Table) <- c(" |","Error Parameters "," |","Error Values","|")
      message("Error Estimation:")
      print(E.Table, right=F, row.names = F)
    }
    message("Weber-Morris Diffusion Model")
    message("Formula: qt = (kWAM*sqrt(t))+C")
    message("First Linear Summary")
    print(summary(lm(firstlineary~firstlinearx)))
    message("Second Linear Summary")
    print(summary(lm(secondlineary~secondlinearx)))
    params1(slp1,int1,slp2,int2)
    predval(x,n.dat)
    error(x)
    xval1 <- seq(min(x),brk.pts[1],length=100)
    xval2 <- seq(brk.pts[1],seg.range[2], length=100)
    xtension1 <-seq(brk.pts[1],1.25*brk.pts[1],length=10)
    xtension2 <-seq(brk.pts[1],0.65*brk.pts[1],length=10)
    eqm1 = slp[1]*xval1+ y.int[1]
    eqm2 = slp[2]*xval2+ y.int[2]
    eqx1 = slp[1]*xtension1+ y.int[1]
    eqx2 = slp[2]*xtension2+ y.int[2]
    plot.dat1 <- data.frame (xval1,eqm1)
    plot.dat2 <- data.frame (xval2,eqm2)
    plot.dat3 <- data.frame (xtension1,eqx1)
    plot.dat4 <- data.frame (xtension2,eqx2)
    theme_set(theme_bw())
    plot <- ggplot(dat, aes(x=x,y=y))+
      geom_line(data=plot.dat1, aes(xval1,eqm1),size=1, color="red")+
      geom_line(data=plot.dat2, aes(xval2,eqm2),size=1, color="blue")+
      geom_line(data=plot.dat3, aes(xtension1,eqx1),size=1, color="red", linetype= "dashed")+
      geom_line(data=plot.dat4, aes(xtension2,eqx2),size=1, color="blue", linetype= "dashed")+
      geom_vline(xintercept = brk.pts,linetype = "dotted",size=1, color="gray60")+
      geom_point(size=2)+
      labs(subtitle="Plot of the qt vs square root of time with piecewise linear Weber-Morris model",
           y="qt",
           x="sqrt(time)",
           title="Weber-Morris Intraparticle Model",
           caption="Created by PUPAK using ggplot2")
    print(plot)
  }
}