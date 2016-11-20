#' Select s0 using the method of Tusher, Tibshirani and Chu (2001)
#' @param beta Effect size estimates
#' @param sds Estimates of sqrt(var(beta))
#'@export
choose_s0 <- function(beta, sds){
  salpha <- as.numeric(quantile(sds, probs=seq(0.01, 0.99, by=0.01)))
  nn <- length(salpha)
  ix <- sapply(sds, FUN=function(w){ sum(w >= c(0, salpha))})

  #v <- as.numeric(by(data=beta/sds, INDICES = ix, FUN = mad, na.rm=TRUE))
  #cvs <- sd(v)/mean(v)
  #tol <- max(min(sds[sds > 0]), 1e-3)

  fct <- function(s0, beta, sds, ix){
    cat(s0, " ")
    xx <- beta/(sds + s0)
    v <- as.numeric(by(data=xx, INDICES = ix, FUN = mad, constant = 1/0.64))
    cv <- sd(v)/mean(v)
    cat(cv, "\n")
    return(cv)
  }

  #W <- optimize(fct, interval=c(log10(min(sds))-1, max(log10(sds))), beta=beta,
  #              sds=sds, ix=ix, tol=0.1)

  zz <- sapply(salpha, FUN=fct, beta=beta, sds=sds, ix=ix)
  return(salpha[which.min(zz)])
}

#'@export
choose_zmin <- function(beta, sds, s0, pos, bandwidth, smoother=c("ksmooth_0", "ksmooth"), zmin_quantile=0.9){
  stopifnot(length(quantile) %in% c(1, 2))
  smoother <- match.arg(smoother)
  if(smoother=="ksmoooth_0"){
    smooth.func <- ksmooth_0
  }else{
    smooth.func <- function(x, y, bandwidth){ ksmooth(x=x, y=y, x.points=x, bandwidth=bandwidth)$y}
  }
  y <- beta/(sds + s0)
  ys <-smooth.func(pos, y, bandwidth)
  if(length(zmin_quantile)==1) return(quantile(abs(ys), zmin_quantile))
  return(quantile(ys, zmin_quantile))
}
