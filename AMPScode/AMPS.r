# -----------------------------------------------------------
# AMPS function in R
# Cornell University
# Daniel Fink 
# Translated into R by Anestis
#
# AMPS constructed with FDR selected 0-order knots and 
# a penalized linear regression spline with knots spaced
# at evenly spaced sample quantiles. Difference operator 
# penalties are applied for all knots. Smoothing parameter
# selection is chosen by minimizing estimated risk.
#
# Dependencies: 
#   VC - Variance Calibration
#   FDRnorm - Sequential FDR
#   quantileknots - Quantile-based knots
#   DifMatrix - Difference operator penalties
#   PLSfit - Penalized Least Squares Fit
# -----------------------------------------------------------

AMPS <- function(x, y, nknots, q, nsim, svect, initsig) {
  n <- length(y)
  
  # Variance Calibration variance estimation
  #-------------------------------------------
  # Initial variance estimation
  # fit <- PsplineDR01(x, y, degree, pknots, [], boundstab) 
  # initsig <- sqrt(fit$sigma2hat)
  VC_result <- VC(x, y, nsim, initsig)
  sigmahat <- VC_result$sigmahat
  vccut <- VC_result$vc_cuty
  mcvc <- VC_result$mcvc
  
  # FDR 0-order knot selection
  #------------------------------
  y1 <- diff(y)
  k1 <- (x[2:n] + x[1:(n - 1)]) / 2
  sigmafdr <- sqrt(2) * sigmahat
  FDR_result <- FDRnorm(y1, k1, sigmafdr, q)
  k1_fdr <- FDR_result$k1_fdr
  mpb <- FDR_result$mpb
  
  # Construct Multi-order Design and Penalty Matrix
  #------------------------------------------------
  zk <- k1_fdr
  Xtemp <- cbind(1, x)  # Add intercept and x as first columns
  lcon <- ncol(Xtemp)
  
  # Zero-order basis for selected knots
  for (alpha in 0) {
    for (i in seq_along(zk)) {
      z <- ((x - zk[i])^alpha) * (x >= zk[i])
      Xtemp <- cbind(Xtemp, z)
    }
  }
  
  # First-order basis with quantile knots
  k <- quantileknots(x, nknots)
  for (alpha in 1) {
    for (i in seq_along(k)) {
      z <- ((x - k[i])^alpha) * (x >= k[i])
      Xtemp <- cbind(Xtemp, z)
    }
  }
  
  # Penalty Matrix
  D1 <- DifMatrix(ncol(Xtemp) - lcon, 1)
  D1 <- t(D1) %*% D1
  D3 <- rbind(
    cbind(diag(0, lcon), matrix(0, lcon, ncol(D1))),
    cbind(matrix(0, nrow(D1), lcon), D1)
  )
  
  # Fit the PLS with X, and D
  # ------------------------------
  fit <- PLSfit(Xtemp, D3, y, sigmahat, svect)
  
  # AMPS Structure
  # ------------------
  ampsfit <- list(
    yhat = fit$yhat,
    estrisk = fit$estrisk,
    U = fit$U,
    z = fit$z,
    bhat = fit$bhat,
    varb = fit$varb,
    X = fit$X,
    D = fit$D,
    y = fit$y,
    sigmahat = fit$sigmahat,
    k1_fdr = k1_fdr,
    mpb = mpb
  )
  
  return(ampsfit)
}

