

choose_s0 <- function(beta, sds){
  salpha <- c(0, as.numeric(quantile(sds, probs=seq(0, 1, by=0.05))))
  nn <- length(salpha)
  ix <- sapply(sds, FUN=function(w){ sum(w >= salpha[-nn])})
  cvs <- c()
  v <- as.numeric(by(data=beta/sds, INDICES = ix, FUN = mad, na.rm=TRUE))
  cvs <- c(sd(v)/mean(v))
  tol <- max(min(sds[sds > 0]), 1e-3)

  fct <- function(s0, beta, sds, ix){
    cat(s0, " ")
    xx <- beta/(sds + s0)
    v <- as.numeric(by(data=xx, INDICES = ix, FUN = mad, constant = 1/0.64))
    cv <- sd(v)/mean(v)
    cat(cv, "\n")
    return(cv)
  }

  W <- optimize(fct, interval=range(sds), beta=beta,
                sds=sds, ix=ix, tol=tol)
  return(W$minimum)
}

dnase_choose_s0_z0 <- function(file.name, pheno.file,
                               maxit=50, z0_quantile=0.9,
                               has.stats=FALSE, tempfile="temp.RData"){

  dat <- read_delim(file.name, delim=" ")
  n <- nrow(dat)
  p <- ncol(dat)
  if(has.stats){
    dat <- dat[,-p]
    p <- p-1
  }


  X <- read_delim(pheno.file, col_names=FALSE, delim=" ")
  x <- as.numeric(X[match(names(dat)[-1], X[,1]),2])

  tt <- huber_helper(dat[, -1], x, maxit=50)
  miss.ix <- which(is.na(tt[2,]))
  tt <- tt[, -miss.ix]

  #save(tt, file=tempfile)
  s= choose_s0(beta=tt[1,],sds=tt[2,])


  #xs <- ksmooth(x=dat[,1], y=tt[1,]/(tt[2,] + s), bandwidth=bandwidth, x.points=dat[,1])$y
  xs <- ksmooth_0(x=dat[-miss.ix,1], y=tt[1,]/(tt[2,] + s), bandwidth=bandwidth)
  save(tt, xs, file=tempfile)
  z0 <- as.numeric(quantile(abs(xs), probs=z0_quantile))
  return(list("sm_stat"=xs, "s0"=s, "z0"=z0))
}
