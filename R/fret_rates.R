
#'@title Calculate lambda and fdr for each peak
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param file.list List of files with output from fret_rates_prelim
#'@param n.perm Number of permutations
#'@param fdr.max Maximum fdr to keep data for
#' @return An object that can be passed to fret_thresholds
#'@export
fret_rates <- function(file.list, fdr.max=0.8, temp.dir=".", write.temp=FALSE){
  #Get info from files
  R <- getobj(file.list[1])
  max1 <- R$m1
  max.perm <- R$mperm
  zmin <- R$zmin
  R$seg.bounds$name <- paste0(R$seg.bounds$chr, ".", 1:nrow(R$seg.bounds))
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
    R$seg.bounds$name <- paste0(R$seg.bounds$chr, ".", 1:nrow(R$seg.bounds))
    segment.bounds <- rbind(segment.bounds, R$seg.bounds)
  }


  #Check inputs
  s <- length(zmin)
  stopifnot(s %in% c(1, 2))

  K <- nrow(segment.bounds)
  cat("There are ", K, " segments total.\n")

  max1 <- max1[order(max1$lambda_perbase, decreasing=FALSE), ]
  nbp <- matrix(rep(segment.bounds$nbp, s), byrow=TRUE, nrow=s)
  mlp_ix <- grep("max_lambda", names(segment.bounds))

  max.lambda <- t(segment.bounds[,mlp_ix]*nbp)

  max1$lambda <- sapply(max1$lambda_perbase, FUN=function(r){
    sum(pmin(r*nbp, max.lambda))
  })
  max1$fdr <- max1$lambda/(1:nrow(max1))

  return(list("max1"=max1, "max.lambda.pb"=max.lambda.pb,
              "max.perm"=max.perm, "nbp"=nbp, "zmin"=zmin))
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

match_segments <- function(chr, pos, segment.bounds,
                           parallel=FALSE, cores=parallel::detectCores()-1){
  chroms <- unique(chr)
  #o <- order(segment.bounds$chr, segment.bounds$start, decreasing=FALSE)
  #stopifnot(all(o==1:nrow(segment.bounds)))
  ix <- rep(NA, length(chr))
  if(parallel){
    cl <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(cl))
  }
  for(c in chroms){
    cat(c, "..")
    ix_dat <- which(chr==c)
    ix_sb <- which(segment.bounds$chr==c)
    stopifnot(all(pos[ix_dat] <= max(segment.bounds$stop[ix_sb])))
    if(parallel){
      slocal <- parSapply(cl, pos[ix_dat], FUN=function(p){
        which.min(c(segment.bounds$start[ix_sb], Inf) <= p)-1
      })
    }else{
      slocal <- sapply(pos[ix_dat], FUN=function(p){
        which.min(c(segment.bounds$start[ix_sb], Inf) <= p)-1
      })
    }
    stopifnot(all(pos[ix_dat]) <= segment.bounds$stop[slocal])
    stopifnot(length(ix[ix_dat])== length(ix_sb[slocal]))
    ix[ix_dat] <- ix_sb[slocal]
  }
  return(ix)
}
