# -----------------------------------------------------------------
# PLS MIN
# Objective Function for Penalized Least Squares 
# see Beran's PLS (2001) paper, equation 2.15 
# 
# lambda = vector of eigenvalues (order)
# z2 = vector of transformed data
# ghat= vector
# Note: make sure that these vectors are in the
# same order.
#
# It is assumed!!! that none of the z2 are zero-
# otherwise ghat is infinity!
#
# a=fminbnd('plsmin',-100,100,optimset('Display','off'),lambda,ghat,z2);
# -----------------------------------------------------------------
	plsmin <- function(nu, lambda, ghat, z2) {
	  # Compute the shrinkage factors
	  f <- 1 / (1 + nu * lambda)
  
	  # Calculate the mean squared error
	  mse <- mean(((f - ghat)^2) * z2)
  
	  return(mse)
	}