

ndiff <- function(n, d = 1) {
# Construct the matrix for n-th differences
  if (d == 1)
    {D <- diff(diag(n))}
  else
    {D <- diff(ndiff(n, d - 1))}
  D}

tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)


bbase <- function(x, knots0, knots1, deg){
# Construct regression-spline basis
    ones<-rep(1,length(x))
    P <- outer(x, knots0, tpower, 0)
    P1 <- outer(x, knots1, tpower, deg)
    if (deg==1) B <- cbind(ones,x,P,P1)
    if (deg==2) B<-cbind(ones,x,x^2,P,P1)
    B }

