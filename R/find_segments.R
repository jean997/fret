

#' Use variance of smoothed  permutation test statistics to define interval endpoints
#' @description Use variance of smoothed  permutation test statistics to define interval endpoints.
#' Intervals are placed so that high variance regions are separated from each other and from low variance
#' regions.
#' @param vv Variance of smoothed permutation test statistics
#' @param pos Position vector. Should be the same length as vv
#' @param min.length Minimum interval length
#' @param z0, z Two thresholds used to define ``high variance regions''.
#' @details Currently, find_segments will smooth vv using a bandwidth of min.length/2.
#' If you don't supply z and z0, z will be the 75th percentlie of these smoothed values and z0 will be
#' the 25th percentile.
#' @return 2 by K matrix listing interval endpoints.
#'@export
find_segments <- function(vv, pos, min.length, z0=NULL, z=NULL,
                          bandwidth=NULL, q=0.05){
  stopifnot(length(vv)==length(pos))
  if(is.null(bandwidth)) bandwidth <- min.length/2
  stopifnot(bandwidth > 0)

  vvs <- ksmooth_0(x=pos, y=vv, xout=pos, bandwidth=bandwidth)
  if(is.null(z0) & is.null(z)){
    z0 <- quantile(vvs, 0.5)
    z <- quantile(vvs, 1-q)
  }

  stopifnot(z0 > 0 & z >= z0)
  stopifnot(all(vvs >= 0))
  if(all(vvs < z)){
    strts <- seq(min(pos), max(pos), by=min.length)
    stps <- strts + min.length-1
    stps[length(stps)] <- max(pos)
    return(cbind(strts, stps))
  }
  ivls <- excursions(vvs, z0)
  #Max stat value inside each excurion
  bpoints <- apply(ivls, MARGIN=1, FUN=function(iv){
    m <- max(vvs[iv[1]:iv[2]])
    if(m >= z) return(iv)
    return(c())})
  bpoints <- unlist(bpoints)
  bp <- matrix(bpoints, ncol=2, byrow=TRUE)
  bp[,1] <- pos[bp[,1]]
  bp[,2] <- pos[bp[,2]]

  strts <- c(1)
  stps <- c()
  ix <- 1
  n <- 1
  while(ix <= nrow(bp)){
    if(bp[ix,1]-strts[n] >= min.length){
      #Low variance region long enough to make up entire segment by itself
      stps <- c(stps, bp[ix, 1]-1)
      strts <- c(strts, bp[ix, 1])
      n <- n+1
    }else if(bp[ix, 2]-strts[n] + 1 >= min.length){
      stps <- c(stps, bp[ix, 2])
      strts <- c(strts, bp[ix, 2] + 1)
      n <- n+1
      ix <- ix + 1
    }else{
      n.needed <- min.length - (bp[ix, 2]-strts[n] + 1 )
      if(n==1) dist.prev=0
        else dist.prev <- min( stps[n-1]-strts[n-1] + 1 - min.length, stps[n-1] - bp[ix-1, 2])
      dist.end <- max(pos)-bp[ix,2]

      if(dist.prev + dist.end < n.needed){
        n <- n-1
        stps[n] <- max(pos)
      }else if(dist.prev >=floor(n.needed/2) & dist.end >= ceiling(n.needed/2)){
        nleft <- floor(n.needed/2)
        nright <- n.needed -nleft
        strts[n] <- strts[n]-nleft
        stps[n-1] <- strts[n]-1
        stps[n] <- bp[ix,2] + nright
      }else if(dist.prev >= floor(n.needed/2)){
        nright <- dist.end
        nleft <- n.needed-nright
        strts[n] <- strts[n]-nleft
        stps[n-1] <- strts[n]-1
        stps <- c(stps, bp[ix,2] + nright)
      }else{
        nleft <- dist.prev
        nright <- n.needed -dist.prev
        strts[n] <- strts[n]-nleft
        stps[n-1] <- strts[n]-1
        stps <- c(stps, bp[ix,2] + nright)
      }
      if(any(bp[,1] <= stps[n] & bp[,2] >= stps[n])){
        ix <- which(bp[,1] <= stps[n] & bp[,2] >= stps[n])
        stps[n] <- bp[ix, 2]
      }
      if(any(bp[,1] >= stps[n])) ix <- min(which(bp[,1] >= stps[n]))
        else ix <- nrow(bp) + 1
      n <- n + 1
      strts[n] <- stps[n-1] + 1
    }
  }
  #strts is length n, stps is length n-1
  if(strts[n] > max(pos)){
    strts <- strts[-n]
    n <- n-1
  }
  #Either stps and strts are both length n or
  #strts is length n, stps is length n-1
  if(strts[n]+min.length-1 <= max(pos)){
    stps[n] <- max(pos)
  }else{
    strts <- strts[-n]
    n <- n-1
    if(length(stps) > n) stps <- stps[1:n]
    stps[n] <- max(pos)
  }
  return(cbind(strts, stps))
}
