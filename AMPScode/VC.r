# ----------------------------------------------
# Variance Calibration with First Differences 
# ----------------------------------------------
# Input
#   nsim = # MC trials
#   x, y = vectors of data
#   initsig = initial estimate of sigma
# Output
#   sigmahat = estimated StDev
#   mcvc = MC-VC Thresholding vector
#   vccut = vc cutoff
# ----------------------------------------------

VC <- function(x, y, nsim, initsig) {
  n <- length(y)
  nvc <- nsim
  ydifi <- numeric(nvc)
  ydif <- numeric(nvc)
  
  # Initial St.Dev. Estimate
  boundstab <- 0
  pknots <- min(100, round(2 * sqrt(length(y))))
  degree <- 1
  
  # Initial Estimate of sigma
  sdadd <- initsig
  p_approx <- NULL
  ymindif <- numeric(nvc)
  wcutsy <- numeric(n - 1)
  
  y1 <- diff(y)  # First Differences
  ysort1 <- sort(abs(y1), index.return = TRUE)
  yindex1 <- ysort1$ix
  
  for (j in 1:nvc) {
    # Adding 1 times original noise works well
    e <- rnorm(n, mean = 0, sd = sdadd)
    obsvar <- mean(e^2)
    yvc <- y + e
    y2 <- diff(yvc)  # First Differences
    ysort2 <- sort(abs(y2), index.return = TRUE)
    yindex2 <- ysort2$ix
    
    dif <- numeric(n - 1)
    dif[1] <- 0
    v1 <- numeric(n - 1)
    v2 <- numeric(n - 1)
    
    for (i in 2:(n - 1)) {
      v1[i] <- mean(y1[yindex1[1:i]]^2)
      v2[i] <- mean(y2[yindex2[1:i]]^2)
      dif[i] <- v2[i] - v1[i]
      
      if (dif[i - 1] < obsvar && dif[i] > obsvar) {
        ydifi[j] <- i
        if (ymindif[j] < ydifi[j]) break
      }
      
         
	    if (v2[i - 1] < obsvar && v2[i] > obsvar) {
        ymindif[j] <- i - 1
      }
    }
    
    # MC-VC Thresholding
    wtcutsy <- numeric(n - 1)
    wtcutsy[yindex1[(ydifi[j] + 1):(n - 1)]] <- 1
    wcutsy <- wcutsy + wtcutsy
  }
  
  # Select the Component cutoff
  vc_cuty <- max(round(mean(ydifi)), round(median(ymindif)))
  
  # Variance Calibration Selected Component #
  sigmahat <- sqrt(sum(y1[yindex1[1:vc_cuty]]^2) / (vc_cuty - 1))
  
  mcvc <- wcutsy
  return(list(sigmahat = sigmahat, vc_cuty = vc_cuty, mcvc = mcvc))
}