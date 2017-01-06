#'@import intervals


#'@title Calculate lambda and fdr for each peak
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param file.list List of files with output from fret_stats
#'@param n.perm Number of permutations
#'@param fdr.max Maximum fdr to keep data for
#' @return An object that can be passed to fret_thresholds
#'@export
fret_rates <- function(file.list, fdr.max=0.8){
  #max1, max.perm, n.perm, zmin, segment.bounds,fdr.max=0.8){
  #Get info from files
  R <- getobj(file.list[1])
  max1 <- R$m1
  max.perm <- R$mperm
  zmin <- R$zmin
  segment.bounds <- R$seg.bounds
  n.perm <- R$n.perm
  file.list <- file.list[-1]
  for(f in file.list){
    cat(f, "\n")
    R <- getobj(f)
    stopifnot(R$zmin==zmin)
    stopifnot(R$n.perm==n.perm)
    max1 <- rbind(max1, R$m1)
    max.perm <- rbind(max.perm, R$mperm)
    segment.bounds <- rbind(segment.bounds, R$seg.bounds)
  }


  #Check inputs
  stopifnot(ncol(segment.bounds)==3)
  stopifnot(names(segment.bounds)==c("chr", "start", "stop"))
  s <- length(zmin)
  stopifnot(s %in% c(1, 2))

  K <- nrow(segment.bounds)
  cat("There are ", K, " segments total.\n")
  nbp <- segment.bounds$stop-segment.bounds$start + 1

  #Segment for each peak in max1
  #max1 has columns mx, chr pos iv1 iv2
  max1$segment <- apply(max1[, c("chr", "pos")], MARGIN=1, FUN=function(x){
    s <- which(segment.bounds$chr==as.character(x[1]) &
                 segment.bounds$start <= as.numeric(x[2]) &
                 segment.bounds$stop >= as.numeric(x[2]) )
    if(length(s)==0) return(NA)
    return(s)
  })
  if(any(is.na(max1$segment))){
    cat("Warning: segment.bounds does not include everything in max1\n")
    max1 <- max1[!is.na(max1$segment),]
  }
  #Segment for each peak in perm.maxes
  max.perm$segment <- apply(max.perm[, c("chr", "pos")], MARGIN=1, FUN=function(x){
    s <- which(segment.bounds$chr==as.character(x[1]) &
                 segment.bounds$start <= as.numeric(x[2]) &
                 segment.bounds$stop >= as.numeric(x[2]) )
    if(length(s)==0) return(NA)
    return(s)
  })
  if(any(is.na(max.perm$segment))){
    cat("Warning: segment.bounds does not include everything in perm.maxes\n")
    max.perm <- max.perm[!is.na(max.perm$segment),]
  }
  max.perm$lambda_perbase <- rep(NA, nrow(max.perm))
  max1$lambda_perbase <- rep(NA, nrow(max1))

  max.lambda.pb <- matrix(nrow=s, ncol=K)
  for(i in 1:K){
    if(i %% 1000 == 1) cat(i, "..")
    m1.ix <- which(max1$segment==i)
    perm.ix <- which(max.perm$segment==i)

    if(length(m1.ix)==0 & length(perm.ix)< 2){
      max.lambda.pb[, i] <- 0
      next
    }
    if(length(perm.ix) > 0){
      m <- max.perm$mx[perm.ix]
      o <- order(m, decreasing=TRUE)
      oinv <- match(1:length(m), o)
      ll <- lamtab(mx=m, zmin=zmin, nbp = nbp[i], n.perm=n.perm)
      max.perm$lambda_perbase[perm.ix] <- ll[,2][oinv]
      if(s==1){
        max.lambda.pb[1, i] <- fret:::get_rate_with_thresh(ll, zmin, np=2)
      }else{
        max.lambda.pb[1, i] <- fret:::get_rate_with_thresh(ll, zmin[1], np=2)
        max.lambda.pb[2, i] <- fret:::get_rate_with_thresh(ll, zmin[2], np=2)
      }
      if(length(m1.ix) > 0){
        #Rates
        rts <- sapply(max1$mx[m1.ix], FUN=function(thresh){
            fret:::get_rate_with_thresh(ll, thresh)
        })
        max1$lambda_perbase[m1.ix] <- rts
      }
      max.perm[perm.ix, ] <- max.perm[perm.ix,][o,]
    }else{
      #perm.ix is empty but max1.ix is not
      #i.e. no permutation peaks above z0 but there are some non perm. peaks above zmin
      #i.e. highly significant but cant estimate significance because all perm peaks are too low
      #probably very rare
      max1$lambda_perbase[m1.ix] <- 0
    }

  }

  max1 <- max1[order(max1$lambda_perbase, decreasing=FALSE), ]
  nbp <- matrix(rep(nbp, s), byrow=TRUE, nrow=s)
  max.lambda <- max.lambda.pb*nbp
  max1$lambda <- sapply(max1$lambda_perbase, FUN=function(r){
    sum(pmin(r*nbp, max.lambda))
  })
  max1$fdr <- max1$lambda/(1:nrow(max1))

  return(list("max1"=max1, "max.lambda.pb"=max.lambda.pb,
              "max.perm"=max.perm, "nbp"=nbp, "zmin"=zmin))
}

#'@title Get thresholds for a target FDR level
#'@param obj Object produced by fret_rates
#'@param target.fdr trget fdr level
#'@param stats.files Two column matrix.
#'The first column is chromosome. The second column is the corresponding stats file.
#' @return An object that can be passed to fret_thresholds
#'@export
fret_thresholds <- function(obj, target.fdr, stats.files){

  stopifnot(ncol(stats.files)==2)
  chrs <- unique(obj$max1$chr)
  stopifnot(all(chrs %in% stats.files[,1]))
  stopifnot(length(target.fdr)==1)
  if(target.fdr < min(obj$max1$fdr)) stop("The smallest possible FDR is",
                                              min(obj$max1$fdr),  ".\n")

  K <- dim(obj$max.lambda.pb)[2]
  s <- dim(obj$max.lambda.pb)[1]

  thresholds <- data.frame(matrix(nrow=K, ncol=5))
  names(thresholds) <- c("num.disc", "thresh.pos", "thresh.neg", "chrom", "file")
  #For each segment record 1) # of discoveries 2) pos threshold 3) neg threshold 4) chromosome 5) file

  tot.disc <- sum(obj$max1$fdr <= target.fdr) ###This is the number of discoveries
  #We want to draw thresholds with lambda = target.fdr*total num discoveries
  lam.target <- target.fdr*tot.disc
  tt <- get_thresh_with_rate(obj$max.perm, obj$max.lambda.pb, obj$nbp,
                               lam.target, obj$zmin)
  thresholds$thresh.pos <- tt[1,]
  if(s==1) thresholds$thresh.neg <- -tt[1,]
    else thresholds$thresh.neg <- tt[2,]

  for(j in 1:K) thresholds$num.disc[j] <- sum(obj$max1$segment[1:tot.disc]==j)
  ix <- which(thresholds$num.disc > 0)
  thresholds$chrom[ix] <- obj$max1$chr[match(ix, obj$max1$segment)]
  thresholds$file[ix] <- stats.files[match(thresholds$chrom[ix], stats.files[,1]), 2]
  discoveries <- get_discoveries(max1=obj$max1[1:tot.disc,], thresholds = thresholds)
  ret <- list("thresholds"=thresholds, "discoveries"=discoveries)
  return(ret)
}


get_thresh_with_rate <- function(max.perm, max.lambda.pb, nbp,
                                 lambda, zmin, np=4){
  s <- length(zmin)
  K <- ncol(nbp)
  stopifnot(all(dim(max.lambda.pb)==dim(nbp)))
  thresh <- matrix(nrow=s, ncol=K)
  zmin.mat <- matrix(rep(zmin, each=K), byrow=TRUE, nrow=s )
  lambda.pb <- lambda/sum(nbp)
  while(any(max.lambda.pb < lambda.pb & max.lambda.pb > 0)){
    ix <- which(max.lambda.pb < lambda.pb)
    thresh[ix] <- zmin.mat[ix]
    lambda <- lambda - sum(max.lambda.pb[ix]*nbp[ix])
    nbp[ix] <- 0
    max.lambda.pb[ix] <- -1
    lambda.pb <- lambda/sum(nbp)
  }
  #Positive
  segs <- (1:K)[nbp[1,] > 0]
  zpos <- sapply(segs, FUN=function(k){
    m <- max.perm[max.perm$segment==k & max.perm$mx > 0, c("mx", "lambda_perbase")]
    get_thresh_with_rate1(m, lambda.pb)
  })
  thresh[1, segs] <- zpos
  if(s==1) return(thresh)

  segs <- (1:K)[nbp[2,] > 0]
  zneg <- sapply(segs, FUN=function(k){
    m <- max.perm[max.perm$segment==k & max.perm$mx < 0, c("mx", "lambda_perbase")]
    m$mx <- -1*m$mx
    m <- m[dim(m)[1]:1, ]
    get_thresh_with_rate1(m, lambda.pb)
  })
  thresh[2, segs] <- -1*zneg
  return(thresh)
}

get_thresh_with_rate1 <- function(ll, rate, np=4){
  if(rate < min(ll[,2])){
    ff <- lm(ll[1:np, 1]~log10(ll[1:np, 2]))
    return(ff$coefficients[2]*log10(rate) + ff$coefficients[1])
  }
  return(approx(y=ll[,1], x=log10(ll[,2]),
                xout=log10(rate), yright=min(ll[,1]))$y)
}


get_rate_with_thresh <- function(ll, thresh, np=4){
  if(sum(ll[,1] < 0) < 2 & thresh < 0) return(0)
  if(sum(ll[,1] > 0) < 2 & thresh > 0) return(0)
  if(thresh > max(ll[,1])){
    ff <- lm(log10(ll[1:np, 2])~ ll[1:np, 1])
    return(10^(ff$coefficients[2]*thresh + ff$coefficients[1]))
  }else if(thresh < min(ll[,1]) & thresh < 0){
    n <- nrow(ll)
    ff <- lm(log10(ll[(n-np+1):n, 2])~ ll[(n-np+1):n, 1])
    return(10^(ff$coefficients[2]*thresh + ff$coefficients[1]))
  }
  sgn <- sign(thresh)
  ix <- which(sign(ll[,1])==sgn)

  return(10^(approx(x=abs(ll[ix,1]), y=log10(ll[ix,2]),
                    xout=abs(thresh),  rule=2:1)$y))
}



get_discoveries <- function(max1, thresholds){
  discoveries <- matrix(nrow=nrow(max1), ncol=3)
  ix <- which(thresholds$num.disc > 0)
  chru <- unique(thresholds$chrom[ix])
  for(c in chru){
    file <- unique(thresholds$file[ix][thresholds$chrom[ix]==c])
    stopifnot(length(file)==1)
    stats <- getobj(file)
    ixc <- which(max1$chr==c)
    discoveries[ixc, 3] <- c
    for(i in ixc){
      strt <- which(stats$sts.smooth$pos==max1$start[i])
      stp <- which(stats$sts.smooth$pos==max1$stop[i])
      D <- stats$sts.smooth[strt:stp,]
      tpos <- thresholds$thresh.pos[max1$segment[i]]
      tneg <- thresholds$thresh.neg[max1$segment[i]]
      excr <- excursions(D$ys, c(tpos, tneg))
      discoveries[i,1] <- D$pos[min(excr[,1])]
      discoveries[i,2] <- D$pos[max(excr[,2])]
    }
  }
  return(discoveries)
}
