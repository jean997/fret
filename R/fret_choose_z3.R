#'@import intervals


#'@title Calculate lambda and fdr for each peak
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param max1 Table giving peak height, position, and chromosome for original data
#'@param max.perm Table giving peak height, position, and chromosome for permuted data
#'@param n.perm Number of permutations
#'@param zmin Minimum Threshold (length 1 or 2)
#'@param segment.bounds Data frame with three columns. First column is chromosome.
#'Columns two and three are segment start and stop
#'@param fdr.max Maximum fdr to keep data for
#'@param target.fdr Target fdr values (vector)
#' @return A list with items z, Robs, and fdr.
#'@export
fret_rates <- function(max1, max.perm, n.perm, zmin, segment.bounds,fdr.max=0.8){

  #Check inputs
  stopifnot(ncol(segment.bounds)==3)
  stopifnot(names(segment.bounds)==c("chrom", "start", "stop"))
  s <- length(zmin)
  stopifnot(s %in% c(1, 2))

  K <- nrow(segment.bounds)
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
    #cat(i, "\n")
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
#'@export
fret_thresholds <- function(obj, target.fdr){

  if(any(target.fdr < min(obj$max1$fdr))) cat("Warining: Some requested FDR levels not possible.\n")
  target.fdr <- target.fdr[target.fdr >= min(obj$max1$fdr)]

  segnames <- paste0("segment", 1:ncol(obj$max.lambda.pb))
  K <- dim(obj$max.lambda.pb)[2]
  s <- dim(obj$max.lambda.pb)[1]
  n <- length(target.fdr)
  Robs <- matrix(nrow=n, ncol=K+2)
  z <- array(dim=c(s, n, K+2))
  dimnames(z) <- list(1:s, target.fdr, c("lambda", "fdr", segnames) )
  if(n==0){
    Robs <- data.frame(Robs)
    names(Robs) <- c("lambda", "fdr", segnames)
    ret <- list("Robs"=Robs, "z"=z)
    return(ret)
  }
  Robs[,2] <- target.fdr
  for(i in 1:s) z[i, , 1] <- target.fdr
  for(i in 1:n){
    ff <- target.fdr[i]
    cat(ff, " ")
    ix.lower <- max(which(obj$max1$fdr <= ff)) ###This is the number of discoveries
    ix.upper <- ix.lower + 1
    lam.target <- ff*ix.lower
    Robs[i, 1] <- lam.target
    for(j in 1:s) z[j, i, 2] <- lam.target
    z[, i, -c(1, 2)] <- fret:::get_thresh_with_rate(obj$max.perm, obj$max.lambda.pb, obj$nbp,
                               lam.target, obj$zmin)

    for(j in 1:K) Robs[i, j+ 2] <- sum(obj$max1$segment[1:ix.lower]==j)
  }
  Robs <- data.frame(Robs)
  names(Robs) <- c("lambda", "fdr", segnames)


  ret <- list("Robs"=Robs, "z"=z)
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
  if(sum(ll[,1] > 0) < 2 & thresh < 0) return(0)
  if(sum(ll[,1] < 0) < 2 & thresh > 0) return(0)
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
