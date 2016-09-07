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
    lamtab <- perm.maxes.list[[i]]
    if(s==1){
                    #threhsold           #lambda per base
      rts <- approx(x=lamtab[,1], y=lamtab[,2],
                    xout=m1, yright = -Inf, yleft = max(lamtab[,2]))$y
    }else{
      ixpos_m1 <- which(m1 > 0)
      ixpos_pm <- which(lamtab[,1] > 0)
      rts <- rep(NA, length(m1))
      rts[ixpos_m1] <- approx(x=lamtab[ixpos_pm,1], y=lamtab[ixpos_pm,2],
                    xout=m1[ixpos_m1], yright = 0, yleft = max(lamtab[ixpos_pm,2]))$y
      ixneg_m1 <- which(m1 < 0)
      ixneg_pm <- which(lamtab[,1] <0)
      rts[ixneg_m1]<- approx(x=lamtab[ixneg_pm,1], y=lamtab[ixneg_pm,2],
                        xout=m1[ixneg_m1], yleft = 0, yright = max(lamtab[ixneg_pm,2]))$y

    }
    m1tab <- rbind(m1tab, cbind(m1, rts, rep(i, length(m1))))
  }
  m1tab <- data.frame(m1tab)
  m1tab$sgn <- rep(1, nrow(m1tab))
  m1tab$sgn[m1tab$m1 < 0] <- 2
  names(m1tab)[3] <- "seg"
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
      if(s==2){
        zz <- get_thresh_sgn(perm.maxes.list, lambda.pb, zmin)
        z <- rbind(z, c(lam.target, ff, zz$zpos))
        zneg <- rbind(zneg, c(lam.target, ff, zz$zneg))
      }else{
        zz <- get_thresh_usgn(perm.maxes.list, lambda.pb, zmin)
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

get_thresh_usgn <- function(perm.maxes.list, lambda.pb, zmin, eps=1e-6){
  z <- sapply(perm.maxes.list, FUN=function(m){
    approx(x=m[,2], y=m[,1], xout=lambda.pb, yleft=Inf, yright=zmin)$y
  })
  return(z-eps)
}

get_thresh_sgn <- function(perm.maxes.list, lambda.pb, zmin, eps=1e-6){
  zpos <- sapply(perm.maxes.list, FUN=function(m){
    ix <- which(m[,1] > 0)
    approx(x=m[ix,2], y=m[ix,1], xout=lambda.pb, yleft=Inf, yright=zmin[1])$y
  })
  zneg <- sapply(perm.maxes.list, FUN=function(m){
    ix <- which(m[,1] < 0)
    approx(x=m[ix,2], y=m[ix,1], xout=lambda.pb, yleft=-Inf, yright=zmin[2])$y
  })
  return(list("zpos"=zpos-eps, "zneg"=zneg-eps))
}

