

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
find_segments <- function(vv, pos, min.length, z0=NULL, z=NULL){
  stopifnot(length(vv)==length(pos))
  vvs <- ksmooth(x=pos, y=vv, bandwidth=min.length/2, x.points=pos)$y
  if(is.null(z0) & is.null(z)){
    z0 <- quantile(vvs, 0.25)
    z <- quantile(vvs, 0.75)
  }

  stopifnot(z0 > 0 & z >= z0)
  stopifnot(all(vvs >= 0))
  if(all(vvs < z)){
    strts <- seq(min(pos), max(pos), by=min.length)
    stps <- strts + min.length-1
    stps[length(stps)] <- max(pos)
    return(cbind(strts, stps))
  }
  q0 <-rle( vvs > z0 )
  p0 <- length(q0$lengths)
  #Excursions at z0
  ivls <- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
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
      n <- n + 1
    }else if(bp[ix, 2]-strts[n] + 1 >= min.length){
      stps <- c(stps, bp[ix, 2])
      strts <- c(strts, bp[ix,2]+1)
      n <- n + 1
      ix <- ix + 1
    }else{
      new.end <- strts[n] + min.length -1
      if(any(bp[,1] < new.end & bp[,2] > new.end)){
        ix <- which(bp[,1] < new.end & bp[,2] > new.end)
        new.end <- bp[ix, 2]
      }
      stps <- c(stps, new.end)
      strts <- c(strts, new.end + 1)
      if(any(bp[,1] > new.end)) ix <- min(which(bp[,1] > new.end))
        else ix <- nrow(bp) + 1
      n <- n + 1
    }
  }
  if(strts[n] > max(pos)){
    strts <- strts[-n]
  }else if(strts[n]+min.length-1 <= max(pos)){
    stps <- c(stps, max(pos))
  }else{
    strts <- strts[-n]
    n <- n-1
    stps[n] <- max(pos)
  }
  return(cbind(strts, stps))
}
