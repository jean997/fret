#' Select s0 using the method of Tusher, Tibshirani and Chu (2001)
#' @param beta Effect size estimates
#' @param se Estimates of sqrt(var(beta))
#'@export
choose_s0 <- function(beta, se){

  if(any(is.na(se))){
    nmiss <- sum(is.na(se))
    cat("Warning: ", nmiss, " positions have missing sd.\n")
    ix <- which(is.na(se))
    beta <- beta[-ix]
    se <- se[-ix]
  }

  salpha <- as.numeric(quantile(se, probs=seq(0.01, 0.99, by=0.01)))
  nn <- length(salpha)
  ix <- sapply(se, FUN=function(w){ sum(w > c(0, salpha))})

  fct <- function(s0, beta, se, ix){
    #cat(s0, " ")
    xx <- beta/(se + s0)
    v <- as.numeric(by(data=xx, INDICES = ix, FUN = mad, constant = 1/0.64))
    #cv <- sd(v)/mean(v)
    #m <- max(abs(v-median(v))/mad(v, constant=1))
    return(v)
  }

  zz <- sapply(salpha, FUN=fct, beta=beta, se=se, ix=ix)
  return(salpha[which.min(zz)])
}

#'@export
choose_zmin <- function(beta, se, s0, pos, bandwidth,
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

  if(any(is.na(se))){
    nmiss <- sum(is.na(se))
    cat("Warning: ", nmiss, " positions have missing sd -- will be treated as zero if s0 > 0.\n")
    ix <- which(is.na(se))
    if(s0 > 0){
      se[ix] <- 0
    }else{
      beta <- beta[-ix]
      se <- se[-ix]
      pos <- pos[-ix]
    }
  }

  y <- beta/(se + s0)
  ys <-smooth.func(pos, y, bandwidth)
  if(length(zmin_quantile)==1) return(quantile(abs(ys), zmin_quantile))
  return(quantile(ys, zmin_quantile))
}
