#'@useDynLib fret
#'@importFrom Rcpp sourceCpp

#'@title Box kernel smoother for integer positions
#'@description Box kernel for use with DNase-seq and similar data. Assumes missing integer positions have observations of 0.
#'@param x Vector of positions (integers)
#'#'@param y Vector of observations to be smoothed
#'@param xout Vector of output positions
#'@param bandwidth Smoother bandwidth. Should be an odd integer
#'@param chunksize Size of chunks to break xout up into. Speeds computation.
#'@return Vector of smoothed values with length equal to length of xout.
#'@export
ksmooth_0_jean <- function(x, y, xout, bandwidth, chunksize){
  stopifnot(length(x)==length(y))
  if(!all(floor(x)==x) | !all(floor(xout)==xout)) stop("ksmooth_0 should only be used with integer positions\n")
  if(!floor(bandwidth)==bandwidth) stop("Please use an odd integer bandwidth with ksmooth_0")
  if(any(is.na(y))) stop("No missing values please.\n")

  if(bandwidth %% 2 == 0){
    cat("Warning: ksmooth_0 must use odd bandwidth. Replacing ", bandwidth,
        " with ", bandwidth + 1, "\n")
    bandwidth <- bandwidth + 1
  }
  if(missing(chunksize)){
    yout <- ksmooth_0_cpp(x, y, xout, bandwidth)
    return(yout)
  }

  yout <- ksmooth_0_stitch(x, y, xout, bandwidth, chunksize)
  return(yout)
}

#'@title Box kernel smoother for integer positions
#'@description Box kernel for use with DNase-seq and similar data. Assumes missing integer positions have observations of 0.
#'@param x Vector of positions (integers)
#'#'@param y Vector of observations to be smoothed
#'@param xout Vector of output positions
#'@param bandwidth Smoother bandwidth. Should be an odd integer
#'@return Vector of smoothed values with length equal to length of xout.
#'@export
ksmooth_0 <- function(x, y, xout, bandwidth){

  stopifnot(length(x)==length(y))
  if(!all(floor(x)==x)) stop("Only integer positions please.\n")
  if(missing(xout)){
    xout <- x
  }
  if(!all(floor(xout)==xout)) stop("Only integer positions please.\n")
  if(!floor(bandwidth)==bandwidth) stop("Please use an odd integer bandwidth with ksmooth_0")
  if(any(is.na(y))) stop("No missing values please.\n")

  if(bandwidth %% 2 == 0){
    warning(paste0("ksmooth_0 must use odd bandwidth. Replacing ", bandwidth,
        " with ", bandwidth + 1, "\n"))
    bandwidth <- bandwidth + 1
  }

  margin <- ceiling(bandwidth/2)
  xlong <- (min(x)-margin):(max(x) + margin)
  ylong <- rep(0, length(xlong))
  ylong[match(x, xlong)] <- y
  if(missing(xout)) xout <- min(x):max(x)
  yout <- ksmooth(xlong, ylong, x.points=xout, bandwidth=bandwidth)$y
  return(yout)
}
