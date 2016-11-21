#' Select s0 using the method of Tusher, Tibshirani and Chu (2001)
#' @param beta Effect size estimates
#' @param sds Estimates of sqrt(var(beta))
#'@export
choose_s0 <- function(beta, sds){
  salpha <- as.numeric(quantile(sds, probs=seq(0.01, 0.99, by=0.01)))
  nn <- length(salpha)
  ix <- sapply(sds, FUN=function(w){ sum(w >= c(0, salpha))})

  if(any(is.na(sds))){
    nmiss <- sum(is.na(sds))
    cat("Warning: ", nmiss, " positions have missing sd.\n")
    ix <- which(is.na(sds))
    beta <- beta[-ix]
    sds <- sds[-ix]
  }

  fct <- function(s0, beta, sds, ix){
    cat(s0, " ")
    xx <- beta/(sds + s0)
    v <- as.numeric(by(data=xx, INDICES = ix, FUN = mad, constant = 1/0.64))
    cv <- sd(v)/mean(v)
    cat(cv, "\n")
    return(cv)
  }

  zz <- sapply(salpha, FUN=fct, beta=beta, sds=sds, ix=ix)
  return(salpha[which.min(zz)])
}

#'@export
choose_zmin <- function(beta, sds, s0, pos, bandwidth,
                        smoother=c("ksmooth_0", "ksmooth"), zmin_quantile=0.9,
                        stitch=1e5, parallel=TRUE){
  stopifnot(length(quantile) %in% c(1, 2))
  smoother <- match.arg(smoother)
  if(smoother=="ksmoooth_0"){
    smooth.func <- function(x, y, bandwidth){
      ksmooth_0(x, y, bandwidth, stitch=stitch, parallel=parallel)
    }
  }else{
    smooth.func <- function(x, y, bandwidth){
      ksmooth(x=x, y=y, x.points=x, bandwidth=bandwidth)$y
    }
  }

  if(any(is.na(sds))){
    nmiss <- sum(is.na(sds))
    cat("Warning: ", nmiss, " positions have missing sd -- will be treated as zero if s0 > 0.\n")
    ix <- which(is.na(sds))
    if(s0 > 0){
      sds[ix] <- 0
    }else{
      beta <- beta[-ix]
      sds <- sds[-ix]
      pos <- pos[-ix]
    }
  }

  y <- beta/(sds + s0)
  ys <-smooth.func(pos, y, bandwidth)
  if(length(zmin_quantile)==1) return(quantile(abs(ys), zmin_quantile))
  return(quantile(ys, zmin_quantile))
}
