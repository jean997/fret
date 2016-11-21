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
ksmooth_0 <- function(x, y, xout, bandwidth, stitch=NULL, parallel=FALSE, cl=NULL){
  if(!all(floor(x)==x) | !all(floor(xout)==xout)) stop("ksmooth_0 should only be used with integer positions\n")
  if(!floor(bandwidth)==bandwidth) stop("Please use integer bandwidth with ksmooth_0")
  if(any(is.na(y))) stop("No missing values please.\n")

  if(bandwidth %% 2 == 0){
    cat("Warning: ksmooth_0 must use odd bandwidth. Replacing ", bandwidth,
        " with ", bandwidth + 1, "\n")
    bandwidth <- bandwidth + 1
  }
  if(is.null(stitch)){
    yout <- ksmooth_0_cpp(x, y, xout, bandwidth)
    return(yout)
  }

  strts2 <- seq(1, length(x), by=stitch)
  strts1 <- c(1, strts2[-1]-bandwidth)
  stps1 <- c(strts2[-1] + bandwidth, length(x))
  stps2 <- c(strts2[-1]-1, length(x))
  N <- length(strts2)

  if(!parallel){
    yout <- unlist(lapply(1:N, FUN=function(ix){
      ksmooth_0_cpp(x[strts1[ix]:stps1[ix]],
                  y[strts1[ix]:stps1[ix]], xout[strts2[ix]:stps2[ix]], bandwidth)
    }))
    return(yout)
  }
  if(is.null(cl)){
    no_cores <- detectCores()-1
    cl <- makeCluster(no_cores, type="FORK")
  }
  yout <- unlist(parLapply(cl, 1:N, FUN=function(ix){
    ksmooth_0_cpp(x[strts1[ix]:stps1[ix]],
                  y[strts1[ix]:stps1[ix]], xout[strts2[ix]:stps2[ix]], bandwidth)
  }))
  stopCluster(cl)
  return(yout)
}

