#'@import intervals


#'@title Find thresholds for a range of lambda values. Calculate FDR.
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param max1.list List of peak heights for each segment.
#'@param perm.maxes.list a list of matrices giving peak heights for permutations. Each matrix has two columns. First column is the threshold level,
#'second column is lambda per bp (num of peaks less than threshold divided by number of permutations and number of base-pairs).
#'@param nbp Number of base-pairs per segment
#'@param zmin
#'@param fdr.max
#'@param seg.names
#'@param target.fdr
#' @return A list with items z, Robs, and fdr.
#'@export
fret_choose_z2 <- function(max1.list, perm.maxes.list, nbp, zmin, fdr.max=0.8,
                            seg.names=NULL, target.fdr=NULL){

  stopifnot(length(perm.maxes.list)==length(nbp))
  stopifnot(length(perm.maxes.list)==length(max1.list))

  s <- length(zmin)

  K <- length(max1.list)
  if(is.null(seg.names)) seg.names <- paste0("segment", 1:K)
  m1tab <- matrix(nrow=0, ncol=3)
  for(i in 1:K){
    m1 <- sort(max1.list[[i]], decreasing=TRUE)
    if(length(m1)==0) next
    lamtab <- perm.maxes.list[[i]]
    rts <- sapply(m1, FUN=function(thresh){
      get_rate_with_thresh(lamtab, thresh)
    })
    cat(i, " ", class(rts), " ")
    m1tab <- rbind(m1tab, cbind(m1, rts, rep(i, length(m1))))
  }
  m1tab <- data.frame(m1tab, row.names=NULL)
  names(m1tab)[3] <- "seg"
  m1tab$sgn <- rep(1, nrow(m1tab))
  m1tab$sgn[m1tab$m1 < 0] <- 2

  m1tab <- m1tab[order(m1tab$rts, decreasing=FALSE), ]
  m1tab$lambda <- m1tab$rts*sum(nbp)*s
  m1tab$fdr <- m1tab$lambda/(1:nrow(m1tab))
  if(all(m1tab$fdr > fdr.max)){
    cat("No excursions with fdr value less than ", fdr.max, "\n")
    ret <- list("Robs"=NULL, "z"=NULL, "zneg" = NULL, "nbp"=nbp, "m1tab"=m1tab)
    return(ret)
  }
  ix <- min(nrow(m1tab),  max(which(m1tab$fdr <= fdr.max)) + 1)
  z <- Robs <- matrix(nrow=ix, ncol=K+2)
  z[,1] <- Robs[,1] <- m1tab$lambda[1:ix]
  z[,2] <- Robs[,2] <- m1tab$fdr[1:ix]
  if(s==1){
    zneg <- NULL
  }else{
    zneg <- z
  }
  robs <- rep(0, K)
  for(i in 1:ix){
    seg <- m1tab$seg[i]
    robs[seg] <- robs[seg] + 1
    Robs[i,-c(1, 2)] <-robs
    #Thresholds
    if(s==2){
      zz <- get_thresh_sgn(perm.maxes.list, m1tab$rts[i], zmin)
      z[i,-c(1, 2)] <- zz$zpos
      zneg[i, -c(1, 2)] <- zz$zneg
    }else{
      zz <- get_thresh_usgn(perm.maxes.list, m1tab$rts[i], zmin)
      z[i, -c(1, 2)] <- zz
    }
  }
  if(!is.null(target.fdr)){
    target.fdr <- target.fdr[target.fdr >= min(m1tab$fdr)]
    for(ff in target.fdr){
      ix.lower <- max(which(m1tab$fdr <= ff))
      ix.upper <- ix.lower + 1
      lam.target <- ff*ix.lower
      lambda.pb <- lam.target/(sum(nbp)*s)
      if(lambda.pb > m1tab$rts[ix.upper]) stop("Something is wrong!\n")
      Robs <- rbind(Robs, Robs[ix.lower,])
      zz <- get_thresh_with_rate(perm.maxes.list, lambda.pb, zmin)
      if(s==2){
        z <- rbind(z, c(lam.target, ff, zz$zpos))
        zneg <- rbind(zneg, c(lam.target, ff, zz$zneg))
      }else{
        z <- rbind(z, c(lam.target, ff, zz))
      }
    }
  }
  Robs <- data.frame(Robs)
  z <- data.frame(z)
  names(Robs) <- names(z) <- c("lambda", "fdr", seg.names)
  if(s==2){
    zneg <- data.frame(zneg)
    names(zneg) <- names(z)
  }
  ret <- list("Robs"=Robs, "z"=z, "zneg" = zneg, "nbp"=nbp, "m1tab"=m1tab)
  return(ret)
}


get_thresh_with_rate <- function(perm.maxes.list, lambda.pb, zmin, eps=1e-6, np=4){
  if(length(zmin)==1){
    z <- sapply(perm.maxes.list, FUN=function(m){
      get_thresh_with_rate1(m, lambda.pb)
    })
    z <- pmax(z-eps, zmin)
    return(z-eps)
  }
  zpos <- sapply(perm.maxes.list, FUN=function(m){
    ix <- which(m[,1] > 0)
    get_thresh_with_rate1(m[ix,], lambda.pb)
  })
  zpos <- pmax(zpos-eps, zmin[1])
  zneg <- sapply(perm.maxes.list, FUN=function(m){
    ix <- which(m[,1] < 0)
    get_thresh_with_rate1(m[ix,], lambda.pb)
  })
  zneg <- pmin(zneg+eps, zmin[2])
  return(list("zpos"=zpos, "zneg"=zneg))
}

get_thresh_with_rate1 <- function(lamtab, rate, np=4){
  if(rate < min(lamtab[,2])){
    ff <- lm(lamtab[1:np, 1]~log10(lamtab[1:np, 2]))
    return(ff$coefficients[2]*log10(rate) + ff$coefficients[1])
  }
  return(approx(y=lamtab[,1], x=log10(lamtab[,2]),
                xout=log10(rate), yright=min(lamtab[,1]))$y)
}


get_rate_with_thresh <- function(lamtab, thresh, np=4){

  if(thresh > max(lamtab[,1])){
    ff <- lm(log10(lamtab[1:np, 2])~ lamtab[1:np, 1])
    return(10^(ff$coefficients[2]*thresh + ff$coefficients[1]))
  }else if(thresh < min(lamtab[,1])){
    n <- nrow(lamtab)
    ff <- lm(log10(lamtab[(n-np+1):n, 2])~ lamtab[(n-np+1):n, 1])
    return(10^(ff$coefficients[2]*thresh + ff$coefficients[1]))
  }
  sgn <- sign(thresh)
  ix <- which(sign(lamtab[,1])==sgn)

  return(10^(approx(x=lamtab[ix,1], y=log10(lamtab[ix,2]), xout=thresh)$y))
}
