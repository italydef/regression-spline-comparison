 DifMatrix <- function(n, d) {
    # Create the matrix delta (n-1) by n
    d1 <- -1 * rep(1, n - 1) # Vector of -1s
    delta <- diag(1, n, n) + cbind(rep(0,n),rbind(diag(d1), 0)) # Add diagonal shifted by 1
		
    delta <- delta[-n, ] # Remove the last row
  
    if (d == 0) {
      D <- diag(1, n, n) # Identity matrix of size n
    } else if (d == 1) {
      D <- delta
    } else {
      D <- DifMatrix(n - 1, d - 1) %*% delta
    }
  
    return(D)
  }
