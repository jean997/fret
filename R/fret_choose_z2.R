
#' Find thresholds for a range of lambda values. Calculate FDR.
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param max1.list List of peak heights for each segment.
#'@param perm.maxes.list a list of matrices giving peak heights for permutations. Each matrix has two columns. First column is the threshold level,
#'second column is lambda per bp (num of peaks less than threshold divided by number of permutations and number of base-pairs).
#'@param nbp Number of base-pairs per segment
#'@param nlam Number of lambda values to consider
#' @return A list with items z, Robs, and fdr.
#'@export
fret_choose_z2 <- function(max1.list, perm.maxes.list, nbp, fdr.max=0.5,
                           signed=FALSE, seg.names=NULL){

  stopifnot(length(perm.maxes.list)==length(nbp))
  stopifnot(length(perm.maxes.list)==length(max1.list))
  K <- length(max1.list)

  m1tab <- matrix(nrow=0, ncol=3)
  for(i in 1:K){
    m1 <- sort(max1.list[[i]], decreasing=TRUE)
    lamtab <- perm.maxes.list[[i]]
    if(!signed){
                    #threhsold           #lambda per base
      rts <- approx(x=lamtab[,1], y=log10(lamtab[,2]),
                    xout=m1, yright = 0, yleft = max(log10(lamtab[,2])))$y
    }else{
      ixpos_m1 <- which(m1 > 0)
      ixpos_pm <- which(lamtab[,1] > 0)
      rts <- rep(NA, length(m1))
      rts[ixpos_m1] <- approx(x=lamtab[ixpos_pm,1], y=log10(lamtab[ixpos_pm,2]),
                    xout=m1[ixpos_m1], yright = 0, yleft = max(log10(lamtab[ixpos_pm,2])))$y
      ixneg_m1 <- which(m1 < 0)
      ixneg_pm <- which(lamtab[,1] <0)
      rts[ixneg_m1]<- approx(x=lamtab[ixneg_pm,1], y=log10(lamtab[ixneg_pm,2]),
                        xout=m1[ixneg_m1], yright = 0, yleft = max(log10(lamtab[ixneg_pm,2])))$y

    }
    m1tab <- rbind(m1tab, cbind(m1, rts, rep(i, length(m1))))
  }
  m1tab <- data.frame(m1tab)
  m1tab$sgn <- rep(1, nrow(m1tab))
  m1tab$sgn[m1tab$m1 < 0] <- 2
  names(m1tab)[3] <- "seg"
  m1tab <- m1tab[order(m1tab$rts, decreasing=FALSE), ]

  z <- Robs <- matrix(nrow=0, ncol=K+1)
  fdr <- c()
  if(!signed){
    zneg <- NULL
    robs <- rep(0, K)
    zz <- rep(Inf, K)
    lam_seg <- rep(0, K)
    lambda <- 0
    i <- 1
    while(max(fdr) <= fdr.max & i <= nrow(m1tab)){
      seg <- m1tab$seg[i]
      #lambda
      #lam_seg[seg] <- 10^(m1tab$rts[i])*nbp[seg]
      #lam_seg[lam_seg > 0] <- 10^(m1tab$rts[i])*nbp[lam_seg > 0]
      lambda <- 10^(m1tab$rts[i])*sum(nbp)
      #Number of discoveries
      robs[seg] <- robs[seg] + 1
      Robs <- rbind(Robs, c(lambda, robs))
      #Thresholds
      zz[seg] <- m1tab$m1[i]
      z <- rbind(z, c(lambda, zz))
      fdr <- c(fdr, lambda/sum(robs))
      i <- i + 1
    }
  }else{
    zneg <- z
    robs <- rep(0, K)
    zz <- rep(Inf, K)
    zzneg <- rep(-Inf, K)
    lam_seg <- matrix(0, nrow=2, ncol=K)
    lambda <- 0
    i <- 1
    while(max(fdr) <= fdr.max & i <= nrow(m1tab)){
      seg <- m1tab$seg[i]
      sgn <- m1tab$sgn[i]
      #lambda
      #lam_seg[sgn, seg] <- 10^(m1tab$rts[i])*nbp[seg]
      #lam_seg[1, colSums(lam_seg) > 0] <- 10^(m1tab$rts[i])*nbp[colSums(lam_seg) > 0]
      #lam_seg[2, colSums(lam_seg) > 0] <- 10^(m1tab$rts[i])*nbp[colSums(lam_seg) > 0]
      #lambda <- sum(lam_seg)
      #lam_seg[sgn, seg] <- 10^(m1tab$rts[i])
      lambda <- 2*10^(m1tab$rts[i])*sum(nbp)
      #Number of discoveries
      robs[seg] <- robs[seg] + 1
      Robs <- rbind(Robs, c(lambda, robs))

      #Thresholds
      if(sgn ==1 ){
        zz[seg] <- m1tab$m1[i]
        z <- rbind(z, c(lambda, zz))
        zneg <- rbind(zneg, c(lambda, zzneg))
      }else{
        zzneg[seg] <- m1tab$m1[i]
        z <- rbind(z, c(lambda, zz))
        zneg <- rbind(zneg, c(lambda, zzneg))
      }
      fdr <- c(fdr, lambda/sum(robs))
      i <- i+1
    }
  }
  ret <- list("Robs"=Robs, "z"=z, "zneg" = zneg, "fdr"=fdr, "nbp"=nbp, "m1tab"=m1tab)
  return(ret)
}
