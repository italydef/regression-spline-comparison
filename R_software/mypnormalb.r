"pnormal" <- function(x, y, knots0, knots1, bdeg, pord, lambda, plot = T, se = T, pen=1, sigma=1)
{

    # Function pnormal: smooths scatterplot data with regression-splines.
    # Input: 
    #   x = abcissae of data
    #   y = binary response
    #   knots0 = zero order knots for regression-splines
    #   knots1 = bdeg order knots for regression-splines
    #   bdeg = degree of regression-splines
    #   pord = order of difference penalty
    #   lambda = smoothness parameter
    #   plot = plot parameter (T of F)
    #   se = plot parameter (T or F)
    
    # Output: an object of class "pspfit" with the following fields
    #   bdeg = degree of B-splines
    #   cv = cross-validation sum of squares
    #   ed.resid = effective degrees of freedom residuals
    #   effdim = effective dimension P-spline model
    #   family = "gaussian" (like glm object)
    #   lambda = smoothing parameter
    #   link = "identity" (like glm object)
    #   muhat = expected values for y (at x)
    #   mse = standard deviation of errors
    #   nseg = number of B-spline segments on domain from xmin to xmax)
    #   pord = order of difference penalty
    #   x = x as input
    #   xgrid = x grid used for plotting curve
    #   xmin = left boundary of B-spline domain
    #   xmax = right boundary of B-spline domain
    #   y = y as input
    #   ygrid = computed curve on x grid
    
    #
    # Side effect: a plot of (x,y) and the estimated curve (if plot = T) with twice se bands (if se=T).

    #
    # Paul Eilers and Brian Marx, 2003 (c)
    #

# Compute the spline basis design matrix

    m <- length(x)
    B <- bbase(x, knots0, knots1, bdeg)   

# Construct penalty stuff
    if (pen ==1) {
         n <- dim(B)[2]
         P <- sqrt(lambda) * ndiff(n, pord)
         nix <- rep(0, n - pord) 
     }
    if (pen ==2) {
         n <- dim(B)[2]
         aux1<-rep(0,bdeg+1)
         aux2<-rep(1,n-bdeg-1)
         D<- diag(c(aux1,aux2))
         P <- sqrt(lambda) * D
         nix <- rep(0, n) 
     }


# Fit
    if(lambda == 0) {
        f <- lsfit(B, y, intercept = F)
    }
    if(lambda > 0) {
        f <- lsfit(rbind(B, P), c(y, nix), intercept = F)
    }
    h <- hat(f$qr)[1:m]
    beta <- as.vector(f$coef)
    mu <- B %*% beta    

# Cross-validation and dispersion
   
   penwt<-10^seq(-12,12,length=100)
   mm <- length(penwt)
   eps <- .Machine$double.eps
   xx<-t(B)%*%B
   try(R <- chol(xx+1000*eps*diag(NROW(xx))))
   # R is not positive definite
   if (all(is.na(R))) try(R <- chol(xx+2000*eps*diag(NROW(xx))))
   if (all(is.na(R))) try(R <- chol(xx+4000*eps*diag(NROW(xx))))
   if (all(is.na(R))) try(R <- chol(xx+8000*eps*diag(NROW(xx))))
 
   Rinvt <- solve(t(R))
   eig.out <- eigen(Rinvt%*%D%*%t(Rinvt),TRUE)
   U <- eig.out$vectors
   C <- eig.out$values
   Z <- B%*%t(Rinvt)%*%U
   Zy <- c(t(Z) %*% y)
   gcv <- rep(0,mm)
   g <- (Zy^2-sigma^2)/Zy^2 # vector
   

   for (i in 1:mm){
      oneld <- (1/(1+penwt[i]*C))
      gcv[i] <- sum((oneld-g)^2*Zy^2)
   }
   imin <- min((1:length(gcv))[(((gcv==min(gcv))+0)!=0)])
   lambda <- penwt[imin]

# Fit again

    if (pen ==2) {
           P <- sqrt(lambda) * D
           f <- lsfit(rbind(B, P), c(y, nix), intercept = F)
     h <- hat(f$qr)[1:m]
    beta <- as.vector(f$coef)
    mu <- B %*% beta    
}

# Compute curve on grid
    u <- seq(min(x), max(x), length = 100)
    Bu <- bbase(u, knots0, knots1, bdeg)
    zu <- Bu %*% as.vector(f$coef)  

# Plot data and fit
    if(plot) {
        plot(x, y, type="n")
        lines(u, zu, col = 1)
        if(se) {
            varf <- diag(Bu %*% solve(t(B) %*% B + t(P) %*% P) %*% t(Bu))
            sef <- s * sqrt(varf)
            upperu <- zu + 2 * sef
            loweru <- zu - 2 * sef
            lines(u, upperu, lty = 2, col = 4)
            lines(u, loweru, lty = 2, col = 4)
        }
    }
# Return list
    pp <- list(x = x, y = y, muhat = mu,   bdeg = bdeg, pord
         = pord, lambda = lambda, xgrid = u, ygrid = zu, cv = gcv, effdim = sum(h
        ), ed.resid = m - sum(h), family = "gaussian", link = "identity",  pcoef=beta)
    class(pp) <- "pspfit"
    pp
}
