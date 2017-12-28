#'@import MASS
#'@import parallel
#'@import LaF


#'@title Calculate Huber test statistics
#'@description Calculate Huber test statistics with variance inflation constant.
#'@param Y matrix (p x n)
#'@param x trait values
#'@param s0 Additional variance
#'@param k Threshold for huber estimator in multiples of scale parameter.
#'@param maxit Maixum iterations to pass to rlm.
#'@return 3 by p matrix giving coefficient estimates, sd estimates and statistic (including s0 adjustment)
#'@export
huber_stats <- function(Y, x, s0 = 0,  k=1.345, maxit=50){
  B <- apply(Y, MARGIN=1, FUN=function(y){
    f <- rlm(y~x, psi=psi.huber, k=k, scale.est="Huber", maxit=maxit)
    coef <- summary(f)$coefficients
    if(nrow(coef)==1) return(c(NA, NA, 0))
    b1 <- coef[2, 1]
    s <- coef[2, 2]
    if(is.na(s) & s0 > 0) return(c(b1, s, b1/s0))
      else if(is.na(s)) return(c(b1, s, 0))
    return(c(b1, s, b1/(s+s0)))
  })
  return(B)
}

#'@title Calculate Huber test statistics using parallel package
#' @description Calculate Huber test statistics with variance infaltion constant using parallel package
#'@param Y matrix (p x n)
#'@param x trait values
#'@param s0 Additional variance
#'@param k Threshold for huber estimator in multiples of scale parameter.
#'@param maxit Maixum iterations to pass to rlm.
#'@return 3 by p matrix giving coefficient estimates, sd estimates and statistic (including s0 adjustment)
#'@export
huber_stats_parallel <- function(Y, x, cores=parallel::detectCores()-1,
                                 s0 = 0,  k=1.345, maxit=50, digits=Inf){

  cl <- makeCluster(cores, type="FORK")

  B <- parApply(cl, Y, MARGIN=1, FUN=function(y){
    y <- round(y, digits=digits)
    f <- rlm(y~x, psi=psi.huber, k=k, scale.est="Huber", maxit=maxit)
    coef <- summary(f)$coefficients
    if(nrow(coef)==1) return(c(NA, NA, 0))
    b1 <- coef[2, 1]
    s <- coef[2, 2]
    if(is.na(s) & s0 > 0) return(c(b1, s, b1/s0))
      else if(is.na(s)) return(c(b1, s, 0))
    return(c(b1, s, b1/(s+s0)))
  })
  stopCluster(cl)
  return(B)
}

#'@title Calculate OLS test statistics using parallel package
#'@description Calculate OLS test statistics using parallel package
#'@param Y matrix (p x n)
#'@param x trait values
#'@param s0 Additional variance
#'@param k Threshold for huber estimator in multiples of scale parameter.
#'@param maxit Maixum iterations to pass to rlm.
#'@return 3 by p matrix giving coefficient estimates, sd estimates and statistic (including s0 adjustment)
#'@export
lm_stats_parallel <- function(Y, x, cores=parallel::detectCores()-1, s0 = 0,  k=1.345, maxit=50){

  cl <- makeCluster(cores, type="FORK")

  B <- parApply(cl, Y, MARGIN=1, FUN=function(y){
    f <- lm(y~x)
    coef <- summary(f)$coefficients
    if(nrow(coef)==1) return(c(NA, NA, 0))
    b1 <- coef[2, 1]
    s <- coef[2, 2]
    if(is.na(s) & s0 > 0) return(c(b1, s, b1/s0))
      else if(is.na(s)) return(c(b1, s, 0))
    return(c(b1, s, b1/(s+s0)))
  })
  stopCluster(cl)
  return(B)
}
