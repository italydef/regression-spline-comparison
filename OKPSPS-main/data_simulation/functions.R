# data_simulation/functions.R

# Simulated data functions
# Library of piecewise/smooth test functions used to generate OKPSPS
# simulation data. f1-f4 are shared with the top-level functions.R (same
# benchmark cases: jump, twin spikes, smooth bump, three jumps); f5-f11
# below are specific to this simulation set.

# f1: single jump discontinuity at x = 0.5 (flat left of the jump, smooth
# oscillation to the right).
f1 <- function(x) {
  y <- ifelse(x < 0.5, -1.5, 0.25 * sin(x ^ 2 * pi ^ 1.5))
  return(y)
}

# f2: two narrow spikes/peaks (near x = 0.3 and x = 0.65) from the
# 1/(eps + (x-c)^2) terms, joined at x = 0.6.
f2 <- function(x) {
  y <-
    ifelse(x < 0.6, 1 / (0.01 + (x - 0.3) ** 2), 1 / (0.015 + (x - 0.65) **
                                                        2))
  return(y)
}

# f3: smooth, symmetric bump (no discontinuities) - a baseline "easy" case.
f3 <- function(x) {
  y <- 100 / (exp(abs(10 * x - 5))) + ((10 * x - 5) ^ 5) / 500
  return(y)
}

# f4: smooth oscillation with three vertical jumps at x = 0.4, 0.6, 0.8.
f4 <- function(x) {
  y <- sin(15 * x) + 0.3 * x ^ 2

  jump1 <- 0.4
  jump2 <- 0.6
  jump3 <- 0.8

  y[x < jump1] <- y[x < jump1] + 2
  y[x >= jump1 & x < jump2] <- y[x >= jump1 & x < jump2] - 2
  y[x >= jump2 & x < jump3] <- y[x >= jump2 & x < jump3] + 1

  return(y)
}


# f5: smooth logistic sigmoid, steep transition centered at x = 0.4.
f5 <- function(x) {
  result <- 90 / (1 + exp(-100 * (x - 0.4)))
  return(result)
}

# f6: smooth oscillation with a decaying exponential bump near x = 0.
f6 <- function(x) {
  result <- sin(10 * pi * x) + 2 * exp(-50 * x ^ 3) + 2
  return(result)
}

# f7: oscillation with two jumps: a vertical shift at x = 0.4 and a flat
# plateau at -2 between x = 0.4 and x = 0.6.
f7 <- function(x) {
  y <- 3 * sin(x ^ 2 * pi ^ 1.5)

  jump1 <- 0.4
  jump2 <- 0.6

  y[x < jump1] <- y[x < jump1] + 2
  y[x >= jump1 & x < jump2] <-  -2


  return(y)
}

# f8: plain smooth quadratic (no discontinuities).
f8 <- function(x){
  y <- x**2
  return(y)
}

# f9: plain smooth linear function (no discontinuities).
f9 <- function(x){
  y <- -x
  return(y)
}

# f10: single jump discontinuity at x = 0 (step from -2 to 2).
f10 <- function(x){
  y <- ifelse(x < 0.0, -2, 2)
  return(y)
}

# f11: single jump discontinuity at x = 0 (step from -1 to 1).
f11 <- function(x){
  y <- ifelse(x < 0.0, -1, 1)
  return(y)
}
