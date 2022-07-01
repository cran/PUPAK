#' @title Non-Linear Adsorption Kinetic Models Summary
#' @description Summarized results of parameter and error values collected from non-linear adsorption kinetic models, namely: Avrami, Elovich, Fractional Power, Pseudo-First-Order, Pseudo-nth-Order, Pseudo-Second-Order, and Richie's Equation.
#' @param t the numerical value for contact time
#' @param qt the numerical value for the amount adsorbed at time t
#' @param qe the numerical value for the amount adsorbed at equilibrium
#' @param n the numerical value for the richie's equation order of reaction
#' @param sort.by the name of the statistical error parameter in which the models are sorted in either increasing or decreasing order. The only accepted arguments are "RMSE" for Relative Mean Square Error, 'MAE' for Mean Absolute Error, 'MSE' for Mean Squared Error, 'RAE' for Relative Absolute Error, 'AIC' for Akaike Information Criterion, 'BIC' for Bayesian Information Criterion, 'R2' for Coefficient of Determination, and 'SE' for Standard Error Estimate. This argument is case-sensitive, and failure to input the correct value will yield a summary of models in alphabetical order.
#' @import nls2
#' @import stats
#' @import Metrics
#' @import utils
#' @return the summarized error and parameter values from non-linear adsorption kinetic models.
#' @examples
#' \donttest{
#' t <- c(15,30,45,60,75,90,105,120)
#' qt <- c(3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
#' qe <- 4.68
#' aksummary.nl(t,qt,qe,n=NULL,"SE")}
#' @author Jeff Ryan S. Magalong
#' @author Joshua Z. DelaCruz
#' @author Jeann M. Bumatay
#' @author Chester C. Deocaris
#' @export
aksummary.nl <- function(t,qt,qe,n,sort.by){
  errors <- new.env()
  parameters <- new.env()
  summary.data <- new.env()
  t  <- t      ;qt <- qt
  if (missing(qe)){
    qe <- max(qt)
  }
  else if(is.null(qe)){
    qe <- max(qt)
  }
  else{qe <- qe}
  dat <- data.frame(t,qt)
  if(missing(n)){
    n <- NULL
  }else{}
  if(missing(sort.by)){
    s<-"sort.by"
  }else if(is.null(sort.by)){
    s<-"sort.by"
  }else if(isFALSE(sort.by)){
    s<-"sort.by"
  }else{
    s<-sort.by
  }
  Avrami.sum.nl <- function(t,qt,qe){
    x  <-t ;y  <-qt
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    fxnavm <- y ~ (qe* (1- ((exp(-k1*(x^n))))))
    grdavm <- data.frame(k1 = c(0,10),
                         n = c(0,1))
    cc<- capture.output(type="message",
                        fit23 <- try(nls2::nls2(fxnavm,
                                                data = dat,
                                                start = grdavm,
                                                algorithm = "port",
                                                control = list(maxiter = 1000)),
                                     silent=TRUE))
    if(is.null(fit23)==TRUE){
      fit23 <- nls2(fxnavm,
                    data = dat,
                    start = grdavm,
                    algorithm = "plinear-random",
                    control = list(maxiter = 1500))
      parsavm <- as.vector(coefficients(fit23))
      pars_k1 <- parsavm[1L]
      k1min <- pars_k1*0.9  ; k1max <- pars_k1*1.1
      grdavm1 <- data.frame(k1=c(k1min,k1max),
                            n = c(0,1))
      fit23 <- nls2(fxnavm,
                    start = grdavm1,
                    algorithm = "brute-force",
                    control=list(maxiter=1000))
    }else{}
    summary.data$Avrami <- summary(fit23)
    errors$rmse.avrami <- rmse(y,predict(fit23))
    errors$mae.avrami  <- mae(y,predict(fit23))
    errors$mse.avrami  <- mse(y,predict(fit23))
    errors$rae.avrami  <- rae(y,predict(fit23))
    errors$PAIC.avrami <- AIC(fit23)
    errors$PBIC.avrami <- BIC(fit23)
    errors$SE.avrami   <- sqrt((sum((y-predict(fit23))^2))/(n.dat-2))
    parsavm1 <- as.vector(coefficients(fit23))
    pars_k1  <- parsavm1[1L]    ;parameters$avrami.k1 <- parsavm1[1L]
    pars_n   <- parsavm1[2L]    ;parameters$avrami.n  <- parsavm1[2L]
  }
  Elovich.sum.nl <- function(t,qt){
    x <- t ;y <- qt
    fxnem <- y ~ ((1/beta)* log(1 + alpha*beta*x))
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    grdem <- data.frame(alpha = c(0,1000),
                        beta = c(0,10))
    cc<- capture.output(type="message",
                        fit24 <- try(nls2::nls2(fxnem,
                                                data = dat,
                                                start = grdem,
                                                algorithm = "port",
                                                control = list(maxiter = 2000)),
                                     silent=TRUE))
    if(is.null(fit24)==TRUE){
      fit24 <- nls2(fxnem,
                    data = dat,
                    start = grdem,
                    algorithm = "plinear-random",
                    control = list(maxiter = 1000))
      parsem <- as.vector(coefficients(fit24))
      pars_alpha <- parsem[1L]; pars_beta <- parsem[2L]; pars_lin <- parsem[3L]

      alphamin <- pars_alpha*0.9;alphamax = pars_alpha*1.1
      betamin  <- (pars_beta/pars_lin)*0.9;betamax =  (pars_beta/pars_lin)*1.1
      grdem1 <- data.frame(alpha <- c(alphamin,alphamax),
                           beta <- c(betamin,betamax))
      fit24 <- nls2(fxnem,
                    start = grdem1,
                    algorithm = "brute-force",
                    control=list(maxiter=1000))
    }else{}
    summary.data$Elovich<-summary(fit24)
    errors$rmse.elovich <- (rmse(y,predict(fit24)))
    errors$mae.elovich  <- (mae(y,predict(fit24)))
    errors$mse.elovich  <- (mse(y,predict(fit24)))
    errors$rae.elovich  <- (rae(y,predict(fit24)))
    errors$PAIC.elovich <- AIC(fit24)
    errors$PBIC.elovich <- BIC(fit24)
    errors$SE.elovich   <- sqrt((sum((y-predict(fit24))^2))/(n.dat-2))
    parsem <- as.vector(coefficients(fit24))
    parameters$elovich.alpha <- parsem[1L]; pars_alpha <- parsem[1L]
    parameters$elovich.beta <- parsem[2L];pars_beta <- parsem[2L]
  }
  FractionalPower.sum.nl<-function(t,qt,qe){
    x <- t ;y <- qt
    fxnfpm <- y ~ qe*(alpha*(x)^beta)
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    grdfpm <- data.frame(alpha = c(0,1000),
                         beta = c(0,10))
    cc<- capture.output(type="message",
                        fit25 <- try(nls2::nls2(fxnfpm,
                                                 data = dat,
                                                 start = grdfpm,
                                                 algorithm = "port",
                                                 control = list(maxiter = 1000)),
                                      silent=TRUE))
    if(is.null(fit25)==TRUE){
      fit25 <- nls2(fxnfpm,
                     data = dat,
                     start = grdfpm,
                     algorithm = "plinear-random",
                     control = list(maxiter = 1000))
      parsfpm <- as.vector(coefficients(fit25))
      pars_alpha <- parsfpm[1L]; pars_beta<- parsfpm[2L]; pars_lin <- parsfpm[3L]

      alphamin <- pars_alpha*1.25*pars_lin
      alphamax <- pars_alpha*0.75*pars_lin
      betamin  <- pars_beta*0.9
      betamax  <- pars_beta*1.1
      grdfpm1 <- data.frame(alpha=c(alphamin,alphamax),
                            beta=c(betamin,betamax))
      fit25 <- nls2(fxnfpm,
                     start = grdfpm1,
                     algorithm = "brute-force",
                     control=list(maxiter=1000))
    }else{}
    summary.data$FractionalPower<-summary(fit25)
    errors$rmse.fracpow <-(rmse(y,predict(fit25)))
    errors$mae.fracpow  <- (mae(y,predict(fit25)))
    errors$mse.fracpow  <- (mse(y,predict(fit25)))
    errors$rae.fracpow  <- (rae(y,predict(fit25)))
    errors$PAIC.fracpow <- AIC(fit25)
    errors$PBIC.fracpow <- BIC(fit25)
    errors$SE.fracpow   <- sqrt((sum((y-predict(fit25))^2))/(n.dat-2))
    parsfpm1 <- as.vector(coefficients(fit25))
    pars_alpha <- parsfpm1[1L];parameters$fracpow.alpha <- parsfpm1[1L]
    pars_beta<- parsfpm1[2L];parameters$fracpow.beta <- parsfpm1[2L]
  }
  PFO.sum.nl <- function(t,qt,qe){
    x   <- t ;y   <- qt
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    fxnpfo <- y ~ (qe* (1- exp(-k1*x)))
    grdpfo <- data.frame(k1 = c(0,10))
    cc<- capture.output(type="message",
                        fit26 <- try(nls2::nls2(fxnpfo,
                                                 data = dat,
                                                 start = grdpfo,
                                                 algorithm = "port",
                                                 control = list(maxiter = 1500)),
                                      silent=TRUE))
    if(is.null(fit26)==TRUE){
      fit26 <- nls2(fxnpfo,
                     data = dat,
                     start = grdpfo,
                     algorithm = "plinear-random",
                     control = list(maxiter = 1000))
      parspfo <- as.vector(coefficients(fit26))
      pars_k1 <- parspfo[1L]
      k1min   <- pars_k1*0.9
      k1max   <- pars_k1*1.1
      grdpfo1 <- data.frame(k1 =c(k1min,k1max))
      fit26 <- nls2(fxnpfo,
                     start = grdpfo1,
                     algorithm = "brute-force",
                     control=list(maxiter=1000))
    }else{}
    summary.data$PFO<-summary(fit26)
    errors$rmse.pfo <- rmse(y,predict(fit26))
    errors$mae.pfo  <- mae(y,predict(fit26))
    errors$mse.pfo  <- mse(y,predict(fit26))
    errors$rae.pfo  <- rae(y,predict(fit26))
    errors$PAIC.pfo <- AIC(fit26)
    errors$PBIC.pfo <- BIC(fit26)
    errors$SE.pfo   <- sqrt((sum((y-predict(fit26))^2))/(n.dat-2))
    parspfo <- as.vector(coefficients(fit26))
    parameters$pfo.k1 <- parspfo[1L]
  }
  PNO.sum.nl <- function(t,qt,qe){
    x <- t ;y <- qt
    fxnpno <- y ~ qe*(1-(1/((1+((n-1)*kn*x*(qe^(n-1))))^(1/(n-1)))))
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    grdpno <- data.frame(kn = c(0,10),
                         n = c(1,3))
    cc<- capture.output(type="message",
                        fit27 <- try(nls2::nls2(fxnpno,
                                                 data = dat,
                                                 start = grdpno,
                                                 algorithm = "port",
                                                 control = list(maxiter = 1000)),
                                      silent=TRUE))
    if(is.null(fit27)==TRUE){
      fit27 <- nls2(fxnpno,
                     data = dat,
                     start = grdpno,
                     algorithm = "plinear-random",
                     control = list(maxiter = 1000))
      pars <- as.vector(coefficients(fit27))
      pars_kn <- pars[1L]; pars_n <- pars[2L]
      knmin <- pars_kn*0.9
      knmax <- pars_kn*1.1
      nmin  <- pars_n*0.9
      nmax  <- pars_n*1.1
      grdpno2 <- data.frame(kn = c(knmin,knmax),
                            n = c(nmin,nmax))
      fit27 <- nls2(fxnpno,
                     start = grdpno2,
                     algorithm = "grid-search",
                     control=list(maxiter=1000))
    }else{}
    summary.data$PNO<-summary(fit27)
    errors$rmse.pno <- rmse(y,predict(fit27))
    errors$mae.pno  <- mae(y,predict(fit27))
    errors$mse.pno  <- mse(y,predict(fit27))
    errors$rae.pno  <- rae(y,predict(fit27))
    errors$PAIC.pno <- AIC(fit27)
    errors$PBIC.pno <- BIC(fit27)
    errors$SE.pno   <- sqrt((sum((y-predict(fit27))^2))/(n.dat-2))
    parspno1 <- as.vector(coefficients(fit27))
    parameters$pno.kn <- parspno1[1L]; pars_kn <- parspno1[1L]
    parameters$pno.n <- parspno1[2L]; pars_n <- parspno1[2L]
  }
  PSO.sum.nl <- function(t,qt,qe){
    x <- t ;y <- qt
    dat <- data.frame(x,y)
    n.dat  <- nrow(na.omit(dat))
    fxnpso <- y ~ (((qe^2)*k2*x)/(1+(qe*k2*x)))
    grdpso <- data.frame(k2 = c(0,100))
    cc<- capture.output(type="message",
                        fit28 <- try(nls2::nls2(fxnpso,
                                                 data = dat,
                                                 start = grdpso,
                                                 algorithm = "port",
                                                 control = list(maxiter = 1000)),
                                      silent=TRUE))
    if(is.null(fit28)==TRUE){
      fit28 <- nls2(fxnpso,
                     data = dat,
                     start = grdpso,
                     algorithm = "plinear-random",
                     control = list(maxiter = 1000))
      parspso <- as.vector(coefficients(fit28))
      pars_k2 <- parspso[1L]
      k2min   <- pars_k2*0.9 ;k2max <- pars_k2*1.1
      grdpso  <- data.frame(k2=c(k2min,k2max))
      fit28 <- nls2(fxnpso,
                     start = grdpso,
                     algorithm = "brute-force",
                     control=list(maxiter=1000))
    }else{}
    summary.data$PSO<-summary(fit28)
    errors$rmse.pso <- rmse(y,predict(fit28))
    errors$mae.pso  <- mae(y,predict(fit28))
    errors$mse.pso  <- mse(y,predict(fit28))
    errors$rae.pso  <- rae(y,predict(fit28))
    errors$PAIC.pso <- AIC(fit28)
    errors$PBIC.pso <- BIC(fit28)
    errors$SE.pso   <- sqrt((sum((y-predict(fit28))^2))/(n.dat-2))
    parspso1 <- as.vector(coefficients(fit28))
    parameters$pso.k2 <- parspso1[1L]
    pars_k2 <- parspso1[1L]
  }
  Richie.sum.nl <- function(t,qt,qe,n){
    x <- t ;y <- qt ;qe <- qe
    dat <- data.frame(x,y)
    n.dat <- nrow(na.omit(dat))
    EQ1 <- function(x,y,qe){
      fxnre <- y ~ (qe* (1- exp(-a*x)))
      grdre <- data.frame(a = c(0,10))
      cc<- capture.output(type="message",
                          fit29 <- try(nls2::nls2(fxnre,
                                                  data = dat,
                                                  start = grdre,
                                                  algorithm = "port",
                                                  control = list(maxiter = 1000)),
                                       silent=TRUE))
      if(is.null(fit29)==TRUE){
        fit29 <- nls2(fxnre,
                      data = dat,
                      start = grdre,
                      algorithm = "plinear-random",
                      control = list(maxiter = 1000))
        parsre <- as.vector(coefficients(fit29))
        pars_a <- parsre[1L]; pars_lin <- parsre[2L]
        amin <- pars_a*0.9 ;amax <- pars_a*1.1
        grdre1 <- data.frame(a=c(amin,amax))
        fit29 <- nls2(fxnre,
                      start = grdre1,
                      algorithm = "brute-force",
                      control=list(maxiter=1000))
      }else{}
      summary.data$Richie<-summary(fit29)
      errors$rmse.richie <- rmse(y,predict(fit29))
      errors$mae.richie  <- mae(y,predict(fit29))
      errors$mse.richie  <- mse(y,predict(fit29))
      errors$rae.richie  <- rae(y,predict(fit29))
      errors$PAIC.richie <- AIC(fit29)
      errors$PBIC.richie <- BIC(fit29)
      errors$SE.richie   <- sqrt((sum((y-predict(fit29))^2))/(n.dat-2))
      parsre1 <- as.vector(coefficients(fit29))
      parameters$richie.a <- parsre1[1L]
      parameters$richie.n <- 1
    }
    EQ2 <- function(x,y,qe){
      fxnre <- y ~ ((a*qe*x)/(1 + (a*x)))
      grdre <- data.frame(a = c(0,10))
      cc<- capture.output(type="message",
                          fit29 <- try(nls2::nls2(fxnre,
                                                  data = dat,
                                                  start = grdre,
                                                  algorithm = "port",
                                                  control = list(maxiter = 1000)),
                                       silent=TRUE))
      if(is.null(fit29)==TRUE){
        fit29 <- nls2(fxnre,
                      data = dat,
                      start = grdre,
                      algorithm = "plinear-random",
                      control = list(maxiter = 1000))
        parsre <- as.vector(coefficients(fit29))
        pars_a <- parsre[1L]; pars_lin <- parsre[2L]
        amin <- pars_a*0.9; amax <- pars_a*1.1
        grdre1 <- data.frame(a=c(amin,amax))
        fit29 <- nls2(fxnre,
                      start = grdre1,
                      algorithm = "brute-force",
                      control=list(maxiter=1000))
      }else{}
      summary.data$Richie<-summary(fit29)
      errors$rmse.richie <- rmse(y,predict(fit29))
      errors$mae.richie  <- mae(y,predict(fit29))
      errors$mse.richie  <- mse(y,predict(fit29))
      errors$rae.richie  <- rae(y,predict(fit29))
      errors$PAIC.richie <- AIC(fit29)
      errors$PBIC.richie <- BIC(fit29)
      errors$SE.richie   <- sqrt((sum((y-predict(fit29))^2))/(n.dat-2))
      parsre1 <- as.vector(coefficients(fit29))
      parameters$richie.a <- parsre1[1L]
      parameters$richie.n <- 2
    }
    EQn <- function(x,y,qe){
      fxnre <- y ~ qe-(qe*(1+((n-1)*a*x))^(1/1-n))
      grdre <- data.frame(a = c(0,100),
                          n = c(1,10))
      cc<- capture.output(type="message",
                          fit29 <- try(nls2::nls2(fxnre,
                                                  data = dat,
                                                  start = grdre,
                                                  algorithm = "port",
                                                  control = list(maxiter=1000)),
                                       silent=TRUE))
      if(is.null(fit29)==TRUE){
        fit29 <- nls2(fxnre,
                      data = dat,
                      start = grdre,
                      algorithm = "plinear-random",
                      control = list(maxiter = 100))
        parsre <- as.vector(coefficients(fit29))
        pars_a <- parsre[1L];pars_n <- parsre[2L]
        amin <- pars_a*0.9 ;amax <- pars_a*1.1
        nmin <- pars_n*0.9 ;nmax <- pars_n*1.1
        grdre1 <- data.frame(a=c(amin,amax),
                             n=c(nmin,nmax))
        fit29 <- nls2(fxnre,
                      start = grdre1,
                      algorithm = "brute-force",
                      control=list(maxiter=1000))
      }else{}
      summary.data$Richie<-summary(fit29)
      errors$rmse.richie <- rmse(y,predict(fit29))
      errors$mae.richie  <- mae(y,predict(fit29))
      errors$mse.richie  <- mse(y,predict(fit29))
      errors$rae.richie  <- rae(y,predict(fit29))
      errors$PAIC.richie <- AIC(fit29)
      errors$PBIC.richie <- BIC(fit29)
      errors$SE.richie   <- sqrt((sum((y-predict(fit29))^2))/(n.dat-2))
      parsre1 <- as.vector(coefficients(fit29))
      parameters$richie.a <- parsre1[1L];pars_a <- parsre1[1L]
      parameters$richie.n <- parsre1[2L];pars_n <- parsre1[2L]
    }
    EQn2<-function(x,y,qe,n){
      n <- n
      fxnre <- y ~ qe-(qe*(1+((n-1)*a*x))^(1/1-n))
      grdre <- data.frame(a = c(0,100))
      cc<- capture.output(type="message",
                          fit29 <- try(nls2::nls2(fxnre,
                                                  data = dat,
                                                  start = grdre,
                                                  algorithm = "port",
                                                  control = list(maxiter=1000)),
                                       silent=TRUE))
      if(is.null(fit29)==TRUE){
        fit29 <- nls2(fxnre,
                      data = dat,
                      start = grdre,
                      algorithm = "plinear-random",
                      control = list(maxiter = 100))

        parsre <- as.vector(coefficients(fit29))
        pars_a <- parsre[2L]
        amin <- pars_a*0.9;amax <- pars_a*1.1
        grdre1 <- data.frame(a=c(amin,amax))
        fit29 <- nls2(fxnre,
                      start = grdre1,
                      algorithm = "brute-force",
                      control=list(maxiter=1000))
      }else{}
      summary.data$Richie<-summary(fit29)
      errors$rmse.richie <- rmse(x,predict(fit29))
      errors$mae.richie  <- mae(x,predict(fit29))
      errors$mse.richie  <- mse(x,predict(fit29))
      errors$rae.richie  <- rae(x,predict(fit29))
      errors$PAIC.richie <- AIC(fit29)
      errors$PBIC.richie <- BIC(fit29)
      errors$SE.richie   <- sqrt((sum((y-predict(fit29))^2))/(n.dat-2))
      parsre1 <- as.vector(coefficients(fit29))
      parameters$richie.a <- parsre1[1L]
      parameters$richie.n <- n
      pars_a <- parsre1[1L]
    }
    if(missing(n)){
      EQn(x,y,qe)}
    else if(is.null(n)){
      EQn(x,y,qe)}
    else if(isFALSE(n)){
      EQn(x,y,qe)}
    else if(n==1){
      EQ1(x,y,qe)}
    else if(n==2){
      EQ2(x,y,qe)}
    else{
      EQn2(x,y,qe,n)}
  }
  Avrami.sum.nl(t,qt,qe)
  Elovich.sum.nl(t,qt)
  FractionalPower.sum.nl(t,qt,qe)
  PFO.sum.nl(t,qt,qe)
  PNO.sum.nl(t,qt,qe)
  PSO.sum.nl(t,qt,qe)
  Richie.sum.nl(t,qt,qe,n)
  Avrami.k1 <- parameters[["avrami.k1"]]
  Avrami.n  <-  parameters[["avrami.n"]]
  Elovich.alpha <- parameters[["elovich.alpha"]]
  Elovich.beta  <- parameters[["elovich.beta"]]
  FractionalPower.alpha <- parameters[["fracpow.alpha"]]
  FractionalPower.beta  <- parameters[["fracpow.beta"]]
  PFO.k1 <- parameters[["pfo.k1"]]
  PSO.k2 <- parameters[["pso.k2"]]
  PNO.kn <- parameters[["pno.kn"]]
  PNO.n  <- parameters[["pno.n"]]
  Richie.alpha <- parameters[["richie.a"]]
  Richie.n <-parameters[["richie.n"]]
  param1.name <- c("k.Avrami","alpha","alpha","k.PFO","k.PNO","k.PSO","alpha")
  param2.name <- c("n","beta","beta","-","n","-","n")
  param1.val  <- as.numeric(c(Avrami.k1,Elovich.alpha,FractionalPower.alpha,PFO.k1,PNO.kn,PSO.k2,Richie.alpha))
  param2.val  <- as.numeric(c(Avrami.n,Elovich.beta,FractionalPower.beta," ",PNO.n," ",Richie.n))
  param1.val  <- round(param1.val,digits=6)
  param2.val  <- round(param2.val,digits=6)
  param2.val[is.na(param2.val)]<-"-"
  Models <- c("Avrami","Elovich","FractionalPower","PFO","PNO","PSO","Richie")
  rmse.val <- as.numeric(c(errors[["rmse.avrami"]],errors[["rmse.elovich"]],errors[["rmse.fracpow"]],errors[["rmse.pfo"]],errors[["rmse.pno"]],errors[["rmse.pso"]],errors[["rmse.richie"]]))
  mae.val  <- as.numeric(c(errors[["mae.avrami"]],errors[["mae.elovich"]],errors[["mae.fracpow"]],errors[["mae.pfo"]],errors[["mae.pno"]],errors[["mae.pso"]],errors[["mae.richie"]]))
  mse.val  <- as.numeric(c(errors[["mse.avrami"]],errors[["mse.elovich"]],errors[["mse.fracpow"]],errors[["mse.pfo"]],errors[["mse.pno"]],errors[["mse.pso"]],errors[["mse.richie"]]))
  rae.val  <- as.numeric(c(errors[["rae.avrami"]],errors[["rae.elovich"]],errors[["rae.fracpow"]],errors[["rae.pfo"]],errors[["rae.pno"]],errors[["rae.pso"]],errors[["rae.richie"]]))
  PAIC.val <- as.numeric(c(errors[["PAIC.avrami"]],errors[["PAIC.elovich"]],errors[["PAIC.fracpow"]],errors[["PAIC.pfo"]],errors[["PAIC.pno"]],errors[["PAIC.pso"]],errors[["PAIC.richie"]]))
  PBIC.val <- as.numeric(c(errors[["PBIC.avrami"]],errors[["PBIC.elovich"]],errors[["PBIC.fracpow"]],errors[["PBIC.pfo"]],errors[["PBIC.pno"]],errors[["PBIC.pso"]],errors[["PBIC.richie"]]))
  SE.val   <- as.numeric(c(errors[["SE.avrami"]],errors[["SE.elovich"]],errors[["SE.fracpow"]],errors[["SE.pfo"]],errors[["SE.pno"]],errors[["SE.pso"]],errors[["SE.richie"]]))
  rmse.val <- round(rmse.val,digits=6)
  mae.val  <- round(mae.val,digits=6)
  mse.val  <- round(mse.val,digits=6)
  rae.val  <- round(rae.val,digits=6)
  PAIC.val <- round(PAIC.val,digits=6)
  PBIC.val <- round(PBIC.val,digits=6)
  SE.val   <- round(SE.val,digits=6)
  Summary.Table <- data.frame(Models,rmse.val,mae.val,mse.val,rae.val,PAIC.val,PBIC.val,SE.val,param1.name,param1.val,param2.name,param2.val)
  colnames(Summary.Table) <- c("Models","RMSE","MAE","MSE","RAE","AIC","BIC","SE","Parameter1","     ","Parameter2","      ")
  if(s=="RMSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RMSE),]
    top.model <- as.character(Summary.Table$Models[which.min(rmse.val)])
    message("Summary of Non-Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Root Mean Square Error (RMSE)")
  }else if(s=="MAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MAE),]
    top.model <- as.character(Summary.Table$Models[which.min(mae.val)])
    message("Summary of Non-Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Mean Absolute Error (MAE)")
  }else if(s=="MSE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$MSE),]
    top.model <- as.character(Summary.Table$Models[which.min(mse.val)])
    message("Summary of Non-Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Mean Squared Error (MSE)")
  }else if(s=="RAE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$RAE),]
    top.model <- as.character(Summary.Table$Models[which.min(rae.val)])
    message("Summary of Non-Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Relative Absolute Error (RAE)")
  }else if(s=="AIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$AIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PAIC.val)])
    message("Summary of Non-Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Akaike Information Criterion (AIC)")
  }else if(s=="BIC"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$BIC),]
    top.model <- as.character(Summary.Table$Models[which.min(PBIC.val)])
    message("Summary of Non-Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Bayesian Information Criterion (BIC)")
  }else if(s=="SE"){
    Sort.Summary.Table <- Summary.Table[order(Summary.Table$SE),]
    top.model <- as.character(Summary.Table$Models[which.min(SE.val)])
    message("Summary of Non-Linear Adsorption Kinetic Models: List of Error Values and Parameters sorted by Standard Error (SE)")
  }else{
    Sort.Summary.Table <- Summary.Table
    message("Summary of Non-Linear Adsorption Kinetic Models: List of Error Values and Parameters (Unsorted)")
    top.model <- as.character("Unsorted")
  }
  Param.Table <- Sort.Summary.Table[c(-2,-3,-4,-5,-6,-7,-8)]
  Sort.Summary.Table <-Sort.Summary.Table[c(-9,-10,-11,-12)]
  print(Sort.Summary.Table, right=FALSE,row.names = FALSE)
  print(Param.Table, right=FALSE,row.names = FALSE )
  if(as.character(top.model)!= "Unsorted"){
    message("Summary of the most fitted model (",top.model,")",sep="")
    print(summary.data[[top.model]])
  }else{}
}
