#'Calculate Huber test statistics
#'@description Calculate two sample Huber
#'@param Y matrix (p x n)
#'@param x trait values
#'@param s0 Additional variance
#'@param k Threshold for huber estimator in multiples of scale parameter.
#'@param maxit Maixum iterations to pass to rlm.
#'@return 3 by p matrix giving coefficient estimates, sd estimates and statistic (including s0 adjustment)
#'@export
huber_stats <- function(Y, x, s0 = 0,  k=1.345, maxit=20){
  if(length(x) > 1) stop("Not currently implemented for multivariate traits.\n")
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- rlm(y~labs, psi=psi.huber, k=k, scale.est="Huber", maxit=maxit)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2]
    return(c(b1, s, b1/(s+s0)))
  })
  return(B)
}
