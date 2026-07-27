varcalib <- function(y, sigmc, m) {
  d      <- diff(y)
  n      <- length(d)
  dvd    <- rep(0, n)
  thresh <- rep(0, m)
  ds     <- sort(abs(d))
  for (i in seq_len(m)) {
    dstar  <- diff(y + sigmc * rnorm(n + 1))
    dstars <- sort(abs(dstar))
    for (c in seq_len(n))
      dvd[c] <- (sum(dstars[1:c]^2) - sum(ds[1:c]^2)) / (2 * c) - sigmc^2
    thresh[i] <- which.max(dvd * (dvd < 0))
  }
  cthresh <- mean(thresh)
  estsig  <- sum(ds[1:cthresh]^2) / (2 * c)
  list(c = cthresh, sigma = sqrt(estsig), dvd = dvd)
}