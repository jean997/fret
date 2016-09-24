
#'@export
choose_s0 <- function(beta, sds){
  salpha <- c(0, as.numeric(quantile(sds, probs=seq(0, 1, by=0.01))))
  nn <- length(salpha)
  ix <- sapply(sds, FUN=function(w){ sum(w >= salpha[-nn])})

  #v <- as.numeric(by(data=beta/sds, INDICES = ix, FUN = mad, na.rm=TRUE))
  #cvs <- sd(v)/mean(v)
  #tol <- max(min(sds[sds > 0]), 1e-3)

  fct <- function(logs0, beta, sds, ix){
    cat(logs0, " ")
    s0 <- 10^(logs0)
    xx <- beta/(sds + s0)
    v <- as.numeric(by(data=xx, INDICES = ix, FUN = mad, constant = 1/0.64))
    cv <- sd(v)/mean(v)
    cat(cv, "\n")
    return(cv)
  }

  W <- optimize(fct, interval=c(log10(min(sds))-1, max(log10(sds))), beta=beta,
                sds=sds, ix=ix, tol=0.1)

  #zz <- sapply(log10(salpha), FUN=fct, beta=beta, sds=sds, ix=ix)
  return(W$minimum)
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
