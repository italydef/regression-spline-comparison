#=============================
# RC: Matlab Ruppert and Carroll  (2000) srs estimate
# ============================
RC<-function(yy,xx){
  ydwt<-wd(yy)
  nlevs<-log(length(yy),2)-1
  vscale <- mad(accessD(ydwt, level = nlevs))
  hatsigma<-varcalib(diff(yy),vscale,30)$sigma # hatsigma=vscale
  
  setVariable(matlab, x = xx)
  setVariable(matlab, y = yy)
  evaluate(matlab,'fit=srslocalauto(x,y); z=fit.yhat(:,2);')
  z<-getVariable(matlab, "z")
  pred=z$z[,1]
  return(pred)
}

#=============================
# BASP: the Bayesian version of Ruppert and Carroll’s local penalty RC estimate
# ============================
BAPS<-function(yy,xx){
  nIter <- 3000        # number of Iterations for MCMC
  burn_in <- 1000     # burn-in for MCMC
  degree <-2           # degree of the p-spline 
  nknots <- 35         #number of knots
  nsubknots <-8        #number of subknotsknots
  setVariable(matlab, x = xx)
  setVariable(matlab, y = yy)
  setVariable(matlab, nIter = nIter)
  setVariable(matlab, burn_in = burn_in)
  setVariable(matlab, degree = degree)
  setVariable(matlab, nknots = nknots)
  setVariable(matlab, nsubknots = nsubknots)
  evaluate(matlab,'fit=bayespspline(x,y,nknots,nsubknots,degree,nIter,burn_in);z=fit.yhat;')
  z1=getVariable(matlab, "z")
  pred=z1$z[,1]
  return(pred)
}

#=========================
# DMS estimate from library(miscF) Denison algorithm's DMS estimate
# Fit a variety of curves by a sequence of piecewise polynomials. 
# The number and location of knots are determined by the Reversible Jump MCMC method.
# see Denison, D. G. T., Mallick, B. K., Smith, A. F. M.(1998) Automatic Bayesian Curve Fitting JRSSB vol. 60, no. 2 333-350
# ========================
DMS<-function(yy,xx){
  results <- curve.polynomial.rjmcmc(yy, xx, lambda=1, l=2, l0=1)
  pred<-rowMeans(results$fitted.save)
 return(pred)
}

#=========================
# BARS estimate from library(BARS) Dimatteo algorithm's BARS
# Fit a variety of curves by a sequence of piecewise polynomials. 
# The number and location of knots are determined by the Reversible Jump MCMC method.
# see  Dimatteo, I., Genovese, C. R., Kass, R. E.(2001) Bayesian Curve-fitting with Free-knot Splines Biometrika vol. 88, no. 4 1055-1071
# ========================
BARS<-function(yy,xx){
  out <- bars(x,yy, fits=TRUE)
  pred<-as.vector(out$postmeans)
  return(pred)
}

#================================
# Hybrid adaptive spline Estimate (Whaba)
#=====================================
HAS<-function(yy,x){
  # Initialisations for HAS fitting
  xx  <- (x-min(x))/(max(x)-min(x))
  n<-length(x)
  # -------------------------------------------------------
  pt <- (rep(1,n)%o%unique(xx))[,-c(length(unique(xx)))]
  bas.nul <- stdz(cbind(1,xx-.5))
  bas.cub <- stdz(cubic(xx,unique(xx)))
  bas.con <- stdz(1*cbind((xx)-pt)>0)
  bas.lin <- stdz((xx-pt)*(xx-pt>0))
  bas.quad <- stdz((xx-pt)^2*(xx-pt>0))
  # ----------------------------------------------------------
  # Create basis function lists
  baslist.has     <- list(bas.nul,bas.cub)
  # ==========================================
  pen.mbas <- 25
  # ===========================================
  # HAS with cubic representers
  has.obj  <- bsml(yy, baslist.has, maxbas=pen.mbas, method="has")
  pred  <- has.obj$fit
  return(pred)
}

#===========================
# Trend filtering (TF)
#===========================
TF<-function(yy,x){
  a = trendfilter(yy,pos=x,ord=3)
  cv = cv.trendfilter(a, k=5, mod="lambda")
  pred<-predict(a,lambda=cv$lambda.min)$fit
  return(pred)
}

#=================================================
# Enhanced adaptive spline estimate (EAPS)
#==================================================
EAPS<-function(yy,x){
  my_data <- data.frame(x1 = x, y = yy)
  initial_gam_model <- gam(y ~ s(x1, bs = "cr", k = 14), data = my_data, )
  result <- fit_gam_optim(
    gam_model = initial_gam_model,
    data = my_data,
    n_knots = 14,
    alpha = 1e-07,
    smoothing_method = "GCV.Cp",
    max_iterations = 10000
  )
  rm(my_data)
  pred<-result$model$fitted.values
  return(pred)
}

#=========================
# Fast Adaptive penalized splines (FASP) of Wiesenfarth et al. 
#========================
FAPS<-function(yy,x){
   n<-length(x)
   kn<- default.knots(x,n/4) 
   kn.var<- default.knots(kn,30)
   yfit<-asp2(yy~f(x,basis="os",degree=3,knots=kn,adap=FALSE),  niter = 20, niter.var = 200)
   #yfit<-asp2(yy~f(x,basis="os",degree=3,knots=kn,adap=TRUE),returnFit=TRUE)
   pred<-yfit$fitted
   return(pred)
 }

#=========================
# Jang and Oh enhanced adaptive smoothing splines (EAS)
#========================
EAS<-function(yy,x){
  pilot <- pilotQV(x, yy)
  pred<-pilot$p_fitted
  return(pred)
}

#==============================================================
# CPR Parsimonious spline knot selection ############
#==============================================================
CPR<-function(yy,x){
  DF  <- data.frame(x = x, yy=yy)
  cp3 <- cp(yy ~ bsplines(x, iknots=default.knots(x,12),bknots = c(0,1.00001)), data = DF)
  pred<-as.vector(predict(cp3$fit, se.fit = FALSE))
  return(pred)
}

#===========================
# SARS estimates
#===========================
SARS<-function(yy,x){
  ydwt<-wd(yy)
  nlevs<-log(length(x),2)-1
  vscale <- mad(accessD(ydwt, level = nlevs))
  hatsigma<-varcalib(diff(yy),vscale,30)$sigma # hatsigma=vscale
  out<-sars(x,yy,sd=hatsigma)
  pred<-out$fit
  return(pred)
}

#===================================
# A-splines (adaptive P splines) ASPL
#===================================
ASPL<-function(yy,x){
  num.knots <- max(5, min(floor(length(unique(x))/4), 30))
  knots <- quantile(unique(x), seq(0, 1, length = num.knots + 2))[-c(1,(num.knots+2))]
  degree <- 3
  pen <- 10^seq(-5, 5, 0.20) # Sequence of penalty parameters to test
  # Fit the adaptive spline model
  aridge_fit <- aspline(x, yy, knots, pen, degree = degree)
  # Select the best model (e.g., using ebic criterion) and predict
  best_knots <- aridge_fit$knots_sel[[which.min(aridge_fit$ebic)]]
  a_fit <- lm(yy ~ bSpline(x, knots = best_knots, degree = degree))
  pred <- predict(a_fit, data.frame(x = x))
  return(pred)
}

#==========================================
# freeknotsplines FREEK
#==========================================
FREEK<-function(yy,x){
  xy.freekt <- fit.search.numknots(x, yy, degree = 3, minknot = 3, maxknot = 9)
  pred<-fitted.freekt(xy.freekt)
  return(pred)
}

#==========================================
# Using AMPS 
#=========================================
AMP<-function(yy,x){
  nknots<-30
  qq<-0.1
  nsim=20
  svect<-c(1)
  #===========================
  # SARS estimates
  #===========================
  ydwt<-wd(yy)
  nlevs<-log(length(x),2)-1
  vscale <- mad(accessD(ydwt, level = nlevs))
  hatsigma<-varcalib(diff(yy),vscale,30)$sigma # hatsigma=vscale
  initsig<-hatsigma
  fit <- AMPS(x, yy, nknots, q=qq, nsim=nsim, svect, initsig)
   pred<-fit$yhat[[1]]
   return(pred)
 }

#============================================
# SCHACE
#============================================
SCHA<-function(yy,x){
  results_CV <- main.SCHACE(
     yy,
     folds = 5,
     Bdf_set = seq(4, 8, 1),
     tuning = "crossvalidation",
     methods = "MSE",
     clambda = "lambda1se" )
   pred<-results_CV$"predicted y"
  return(pred)
}
