

#' Use variance of smoothed  permutation test statistics to define interval endpoints
#' @description Use variance of smoothed  permutation test statistics to define interval endpoints.
#' Intervals are placed so that high variance regions are separated from each other and from low variance
#' regions.
#' @param vv Variance of smoothed permutation test statistics
#' @param pos Position vector. Should be the same length as vv
#' @param min.length Minimum interval length
#' @param q Quantile for upper limit for definition of high variance region. Currently 0.05
#' @return 2 by K matrix listing interval endpoints.
#'@export
find_segments2 <- function(vv, pos, min.length, q=0.05){
  stopifnot(length(vv)==length(pos))
  stopifnot(all(vv >= 0))

  z0 <- quantile(vv, 0.5)
  z <- quantile(vv, 1-q)

  p <- length(pos)
  stopifnot(length(vv)==p )
  ivls <- excursions(vv, z0)
  #Max stat value inside each excurion
  bpoints <- apply(ivls, MARGIN=1, FUN=function(iv){
    m <- max(vv[iv[1]:iv[2]])
    if(m >= z) return(iv)
    return(c())})
  bpoints <- unlist(bpoints)
  bp <- matrix(bpoints, ncol=2, byrow=TRUE)
  bp[,1] <- pos[bp[,1]]
  bp[,2] <- pos[bp[,2]]

  strts <- c(pos[1])
  stps <- c()
  i <- 1
  n <- 1
  while(max(stps) <= pos[p]){

    new.end <- strts[i] + min.length -1

    if(any((bp[,1] <= new.end & bp[,2] >= new.end))){
      ix <- which(bp[,1] <= new.end & bp[,2] >= new.end)
      stps[i] <- bp[ix,2] + 1
    }else{
      stps[i] <- new.end
    }
    i <- i + 1
    strts[i] <- stps[i-1] + 1
  }
  strts <- strts[-i]
  i <- i-1
  stps[i] <- pos[p]
  if(strts[i]-stps[i]+1 < min.length){
    strts <- strts[-i]
    stps <- stps[-i]
    i <- i-1
    stps[i] <- pos[p]
  }
  return(cbind(strts, stps))
}
