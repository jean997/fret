#'@import MASS
#'@import parallel
#'@import LaF
#'@import stringr
#'@import dplyr
#'@import stringr

huber_stat <- function(y, x, s0){
  suppressWarnings(f <- rlm(y~x, psi=psi.huber, scale.est="Huber", maxit=50))
  coef <- summary(f)$coefficients
  if(nrow(coef)==1) return(c(NA, NA, 0))
  b1 <- coef[2, 1]
  s <- coef[2, 2]
  if(is.na(s) & s0 > 0) return(c(b1, s, b1/s0))
  else if(is.na(s)) return(c(b1, s, 0))
  return(c(b1, s, b1/(s+s0)))
}

huber_resids <- function(y, X, covariates){
  form <- as.formula(paste0( "y ~", paste0(covariates, collapse="+")))
  suppressWarnings(f <- rlm(form, psi=psi.huber, scale.est="Huber", maxit=50, data=X))
  return(f$residuals)
}

lm_stat <- function(y, x, s0){
  f <- lm(y~x)
  coef <- summary(f)$coefficients
  if(nrow(coef)==1) return(c(NA, NA, 0))
  b1 <- coef[2, 1]
  s <- coef[2, 2]
  if(is.na(s) & s0 > 0) return(c(b1, s, b1/s0))
  else if(is.na(s)) return(c(b1, s, 0))
  return(c(b1, s, b1/(s+s0)))
}

lm_resids <- function(y, X, covariates){
  form <- as.formula(paste0("y ~", paste0(covariates, collapse="+")))
  f <- lm(form, data=X)
  return(f$residuals)
}

qp_stat_binary <- function(y, x, s0){
  zero_val <- 1e-11
  c0 <- max(zero_val, sum(y[x==0]))
  c1 <- max(zero_val, sum(y[x==1]))
  n0 <- sum(x==0)
  n1 <- sum(x==1)
  b1 <- log(c1/n1)-log(c0/n0)
  mu = rep(c0/n0, n0+n1)
  mu[x==1] = c1/n1
  phi = 1/(n0 + n1 - 2)* sum( (y-mu)^2/mu)
  s = sqrt(phi)*sqrt(1/c0 + 1/c1)
  return(c(b1, s, b1/(s+s0)))
}

qp_stat_continuous <- function(y, x, s0){
  f <- glm(y~x, family="quasipoisson")
  coef <- summary(f)$coefficients
  if(nrow(coef)==1) return(c(NA, NA, 0))
  b1 <- coef[2, 1]
  s <- coef[2, 2]
  if(is.na(s) & s0 > 0) return(c(b1, s, b1/s0))
  else if(is.na(s)) return(c(b1, s, 0))
  return(c(b1, s, b1/(s+s0)))
}

qp_resids <- function(y, X, covariates){
  form <- as.formula(paste0( "y ~", paste0(covariates, collapse="+")))
  f <- glm(form, family="quasipoisson")
  return(f$residuals)
}

