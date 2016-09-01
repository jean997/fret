
#'Kernel smoother assuming missing positions are 0
#'@description Kernel smoother assuming missing positions are 0
#'@param x Vector of positions
#'@param y Vector of test statistics
#'@param bandwidth Smoother bandwidth
#'@param maxit Maixum iterations to pass to rlm.
#'@return 2 by p matrix. Top row is coefficient estimate. Bottom row is sd estimates.
#'@export
ksmooth_0 <- function(x, y, bandwidth){
  n <- length(y)
  y.out <- sapply(x, FUN=function(xx){
    sum(y[ x <= (xx+bandwidth/2) & x >= (xx - bandwidth/2)])/(bandwidth+1)
  })
  return(y.out)
}
