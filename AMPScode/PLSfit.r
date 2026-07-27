#-------------------------------------------------------
# PLSfit(X,D,sigmahat,svect)
# Input 
#   X=DesignMatrix,,varest
#   D=PenaltyMatrix
#   sigmahat=estimated st.dev. 
#   svect= indicator for shrinkage vector estimator
# Output   
#   yhat = vector of estimates
#   U = orthogonal basis
#   z = canonical response
#   estrisk = the estimated risk
#   bhat = parameter estimates
#   varb = variance of parameter estimates
#-------------------------------------------------------

PLSfit <- function(X, D, y, sigmahat, svect) {
  # Initialize variables
  sigma2hat <- sigmahat^2
  p <- ncol(X)
  
  # Decomposition
  Rb <- t(X) %*% X
  A <- chol(Rb + 10000 * .Machine$double.eps * diag(nrow(Rb)))
  A <- solve(t(A))  # Equivalent to A = inv(A') in MATLAB
  B <- A %*% D %*% A
  svd_result <- svd(B)
  Gamma <- svd_result$u
  lambda <- svd_result$d
  Upls <- X %*% t(A) %*% Gamma
  zpls <- t(Upls) %*% y
  
  # Shrinkage vector estimation
  f <- list()
  estrisk <- numeric(length(svect))
  yhat <- list()
  
  for (nsvect in seq_along(svect)) {
    if (svect[nsvect] == 1) {
      # Estimated Risk Minimization
      z2 <- zpls^2
      ghat <- 1 - sigma2hat / z2
      # Use optimize in R for fminbnd
      nu <- optimize(plsmin, c(0, 10000), lambda = lambda, ghat = ghat, z2 = z2)$minimum
      f[[nsvect]] <- 1 / (1 + nu * lambda)
      
    } else if (svect[nsvect] == 2) {
      # Minimax Bound
      z2 <- zpls^2
      ghat <- 1 - sigma2hat / z2
      nu <- optimize(plsmin, c(0, 10000), lambda = lambda, ghat = ghat, z2 = z2)$minimum
      ftemp <- 1 / (1 + nu * lambda)
      
      rhohat <- (zpls^2) / sigma2hat
      fMM <- pmax(0, (rhohat - 1) / (rhohat + 1))
      f[[nsvect]] <- pmin(ftemp, fMM)
      
    } else {
      stop("Wrong Code for Shrinkage Vector Estimation")
    }
    
    # Fits and Estimated Risk
    yhat[[nsvect]] <- Upls %*% diag(f[[nsvect]]) %*% zpls
    estrisk[nsvect] <- mean(sigma2hat * f[[nsvect]]^2 + 
                              (zpls^2 - sigma2hat) * (1 - f[[nsvect]])^2)
  }
  
  # Compute summaries
  btemp <- t(A) %*% A %*% t(X)
  bUtemp <- btemp %*% Upls
  Bhat <- btemp %*% do.call(cbind, yhat)
  S <- Upls %*% diag(f[[1]]) %*% t(Upls)  # Smoother matrix for the first svect
  n <- length(y)
  sig2 <- sum((y - yhat[[1]])^2) / (n - sum(diag(S)))
  Bhatvc <- sig2 * bUtemp %*% diag(f[[1]]^2) %*% t(bUtemp)
  
  # Return results as a list
  fit <- list(
    yhat = yhat,
    estrisk = estrisk,
    U = Upls,
    z = zpls,
    bhat = Bhat,
    varb = Bhatvc,
    X = X,
    D = D,
    y = y,
    sigmahat = sigmahat
  )
  
  return(fit)
}

