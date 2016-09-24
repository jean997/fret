#'@useDynLib fret
#'@importFrom Rcpp sourceCpp

#'@title Kernel smoother assuming missing positions are 0
#'@description Kernel smoother assuming missing positions are 0
#'@param x Vector of positions
#'@param y Vector of test statistics
#'@param bandwidth Smoother bandwidth
#'@param maxit Maixum iterations to pass to rlm.
#'@return 2 by p matrix. Top row is coefficient estimate. Bottom row is sd estimates.
#'@export
ksmooth_0_old <- function(x, y, xout, bandwidth){
  if(!all(floor(x)==x) | !all(floor(xout)==xout)) stop("ksmooth_0 should only be used with integer positions\n")
  if(!floor(bandwidth)==bandwidth) stop("Please use integer bandwidth with ksmooth_0")
  if(any(is.na(y))) stop("No missing values please.\n")

  if(bandwidth %% 2 == 0){
    cat("Warning: ksmooth_0 must use odd bandwidth. Replacing ", bandwidth,
          " with ", bandwidth + 1, "\n")
    bandwidth <- bandwidth + 1
  }

  n <- length(y)
  y.out <- sapply(xout, FUN=function(xx){
    sum(y[ x <= (xx+bandwidth/2) & x >= (xx - bandwidth/2)])/(bandwidth)
  })
  return(y.out)
}

#'@title Kernel smoother assuming missing positions are 0
#'@description Kernel smoother assuming missing positions are 0
#'@param x Vector of positions
#'@param y Vector of test statistics
#'@param bandwidth Smoother bandwidth
#'@param maxit Maixum iterations to pass to rlm.
#'@return 2 by p matrix. Top row is coefficient estimate. Bottom row is sd estimates.
#'@export
ksmooth_0 <- function(x, y, xout, bandwidth){
  if(!all(floor(x)==x) | !all(floor(xout)==xout)) stop("ksmooth_0 should only be used with integer positions\n")
  if(!floor(bandwidth)==bandwidth) stop("Please use integer bandwidth with ksmooth_0")
  if(any(is.na(y))) stop("No missing values please.\n")

  if(bandwidth %% 2 == 0){
    cat("Warning: ksmooth_0 must use odd bandwidth. Replacing ", bandwidth,
        " with ", bandwidth + 1, "\n")
    bandwidth <- bandwidth + 1
  }
  yout <- rep(0, length(xout))
  fret:::ksmooth_0_cpp(x, y, bandwidth, xout, yout)
  return(yout)
}

