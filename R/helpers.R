small_to_big_z <- function(z, indices) {
  # for internal use only
  
  d <- length(indices)
  if(d <= 1) {
    return(NULL)
  }

  big_z <- matrix(0, nrow = nrow(z), ncol = d*(d - 1) / 2)
  
  col <- 1
  for(i in 1:(d-1)) {
    for(j in (i+1):d) {
      big_z[, col] <- z[, indices[i]] - z[, indices[j]]
      col <- col + 1
    }
  }
  
  big_z
}