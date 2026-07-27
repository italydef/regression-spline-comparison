# ---------------------------------------------------------------------
# False Discovery Rate (FDR)
# Sequential implementation based on Benjamini and Hochberg (1995)
# Translated to R by [Anestis
#
# Input:
#   y1 - vector of responses
#   k - vector of predictors
#   sigmafdr - standard deviation
#   q - FDR rate
# Output:
#   k_fdr - vector of knots
#   mbp - reordered vector
# ---------------------------------------------------------------------

FDRnorm <- function(y1, k, sigmafdr, q) {
  nf <- length(y1)
  if (length(k) != nf) {
    stop("length(k) != length(y1)")
  }
  
  k1 <- k
  lcon <- 0
  
  # Standardize and sort y1
  mb <- sort(abs(y1) / sigmafdr, decreasing = TRUE, index.return = TRUE)
  mbp <- rev(mb$x)  # Reorder to ascending
  mindex <- rev(mb$ix)
  
  # Threshold calculation
  j <- 1
  nc <- 0
  t <- numeric(nf)
  t[1] <- qnorm(1 - q * 1 / nf, mean = 0, sd = 1)
  
  while (t[j] < mbp[j]) {
    nc <- j
    j <- j + 1
    t[j] <- qnorm(1 - q * j / nf, mean = 0, sd = 1)
  }
  
  # Find contiguous sequences of nonzero FDR coefficients
  if (nc > 0) {
    kseries <- sort(mindex[(nf - nc + 1):nf])
    kall <- k[kseries]
    ball <- y1[kseries]
    kl <- length(kseries)
    kdiff <- diff(kseries)
    closk <- which(kdiff != 1)
    openk <- c(1, closk + 1)
    closk <- c(closk, kl)
  } else {
    kall <- ball <- NULL
    openk <- closk <- integer(0)
  }
  
  # Select Fixed Knots based on FDR Coefficients
  aBhat0 <- abs(y1)
  Bhat0 <- y1
  Bhatvc <- diag(1, nf)  # Identity matrix for variance-covariance
  
  k_fdr <- numeric(0)
  b_fdr <- numeric(0)
  nj <- 1
  
  for (j in seq_along(openk)) {
    # Check length of sequence
    sl <- kseries[closk[j]] - kseries[openk[j]] + 1
    if (sl > 1) {
      cortemp <- Bhatvc[(lcon + openk[j]):(lcon + closk[j]), 
                        (lcon + openk[j]):(lcon + closk[j])]
      
      # Compute correlations from covariances
      cortemp2 <- sweep(cortemp, 2, sqrt(diag(cortemp)), FUN = "/")
      cortemp2 <- sweep(cortemp2, 1, sqrt(diag(cortemp)), FUN = "/")
      
      # Take maximum of sequence
      if (mean(diag(cortemp2, 1)) > 0.6) {
        max_idx <- which.max(aBhat0[kseries[openk[j]]:kseries[closk[j]]])
        k_fdr[nj] <- k[kseries[openk[j]] + max_idx - 1]
        b_fdr[nj] <- Bhat0[kseries[openk[j]] + max_idx - 1]
        nj <- nj + 1
      } else {
        # Take the whole sequence
        for (jj in seq_len(sl)) {
          k_fdr[nj] <- k[kseries[openk[j]] + jj - 1]
          b_fdr[nj] <- Bhat0[kseries[openk[j]] + jj - 1]
          nj <- nj + 1
        }
      }
    } else {
      # Take the single FDR knot
      k_fdr[nj] <- k[kseries[openk[j]]]
      b_fdr[nj] <- Bhat0[kseries[openk[j]]]
      nj <- nj + 1
    }
  }
  
  return(list(k_fdr = k_fdr, mbp = mbp))
}
