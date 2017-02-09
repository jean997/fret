#'@title Get thresholds for a target FDR level
#'@param obj Object produced by fret_rates
#'@param target.fdr trget fdr level
#'@param stats.files Two column matrix.
#'The first column is chromosome. The second column is the corresponding stats file.
#' @return An object that can be passed to fret_thresholds
#'@export
fret_thresholds <- function(obj, target.fdr, stats.files){

  chrs <- unique(obj$max1$chr)
  if(!is.null(stats.files)){
    stopifnot(ncol(stats.files)==2)
    stopifnot(all(chrs %in% stats.files[,1]))
  }
  stopifnot(length(target.fdr)==1)
  if(target.fdr < min(obj$max1$fdr)) stop("The smallest possible FDR is",
                                          min(obj$max1$fdr),  ".\n")

  K <- nrow(obj$segment.bounds)
  s <- length(obj$zmin)

  thresholds <- data.frame(matrix(nrow=K, ncol=6))
  names(thresholds) <- c("num.disc", "thresh.pos", "thresh.neg", "chrom", "file", "name")
  #For each segment record 1) # of discoveries 2) pos threshold 3) neg threshold 4) chromosome 5) file

  tot.disc <- max(which(obj$max1$fdr <= target.fdr)) ###This is the number of discoveries
  #We want to draw thresholds with lambda = target.fdr*total num discoveries
  lam.target <- target.fdr*tot.disc
  stopifnot(sum(obj$max1$lambda <= lam.target)==tot.disc)
  segs <- unique(obj$max1$name[1:tot.disc])
  tt <- fret:::get_thresh_with_rate(obj$max.perm, obj$segment.bounds,
                                    lam.target, obj$zmin, np=10, segs=segs)
  thresholds$thresh.pos <- tt[,1]
  if(s==1) thresholds$thresh.neg <- -tt[,1]
  else thresholds$thresh.neg <- tt[,2]

  for(j in 1:K) thresholds$num.disc[j] <- sum(obj$max1$name[1:tot.disc]==obj$segment.bounds$name[j])

  thresholds$name <- obj$segment.bounds$name
  ix <- which(thresholds$num.disc > 0)
  thresholds$chrom[ix] <- obj$max1$chr[match(thresholds$name[ix], obj$max1$name)]
  if(is.null(stats.files)) return(thresholds)
  fl <- tempfile(pattern=paste0("thresh", target.fdr), tmpdir=".")
  cat("temporary file: ", fl, "\n")
  save(thresholds, file=fl)
  cat("Thresholds found. Retrieving discoveries.\n")
  thresholds$file[ix] <- stats.files[match(thresholds$chrom[ix], stats.files[,1]), 2]
  discoveries <- get_discoveries(max1=obj$max1[1:tot.disc,], thresholds = thresholds)
  ret <- list("thresholds"=thresholds, "discoveries"=discoveries)
  return(ret)
}

get_thresh_with_rate <- function(max.perm, segment.bounds,
                                 lambda, zmin, np=10, segs=NULL){
  s <- length(zmin)
  K <- nrow(segment.bounds)
  stopifnot("name" %in% names(segment.bounds))
  stopifnot("nbp" %in% names(segment.bounds))
  mlp_ix <- grep("max_lambda_perbase", names(segment.bounds))
  stopifnot(length(mlp_ix) > 0)

  thresh <- matrix(nrow=K, ncol=s)
  zmin.mat <- t(matrix(rep(zmin, each=K), byrow=TRUE, nrow=s ))
  max.lambda.pb <- as.matrix(segment.bounds[,mlp_ix, drop=FALSE])
  if(s==2){
    nbp <- cbind(segment.bounds$nbp, segment.bounds$nbp)
  }else{
    nbp <- matrix(segment.bounds$nbp, nrow=K)
  }
  lambda.pb <- lambda/sum(nbp)
  #Some segments have a maximum rate of false discoveries that is
  # less than the target rate. In these segments we can set the
  # threshold to zmin and allow a little more false discoveries
  # in all the other segments.
  # Here we remove the segments with low max.lambda.pb
  # and adjust the total lambda accordingly
  while(any(max.lambda.pb < lambda.pb)){
    ix <- which(max.lambda.pb < lambda.pb)
    thresh[ix] <- zmin.mat[ix]
    lambda <- lambda - sum(max.lambda.pb[ix]*nbp[ix])
    nbp[ix] <- 0
    max.lambda.pb[ix] <- Inf
    lambda.pb <- lambda/sum(nbp)
  }
  #Positive/All if s==1
  ix <- which(nbp[,1] > 0)
  keep.segs <- segment.bounds$name[ix]
  if(!is.null(segs)) keep.segs <- intersect(keep.segs, segs)
  zpos <- sapply(keep.segs, FUN=function(k){
    m <- max.perm[max.perm$name==k & max.perm$mx > 0, c("mx", "lambda_perbase")]
    fret:::get_thresh_with_rate1(m, lambda.pb, np=np)
  })
  thresh[match(keep.segs, segment.bounds$name), 1] <- zpos
  if(s==1) return(thresh)
  #Negative (s==2)
  ix <- which(nbp[,2] > 0)
  keep.segs <- segment.bounds$name[ix]
  if(!is.null(segs)) keep.segs <- intersect(keep.segs, segs)
  zneg <- sapply(keep.segs, FUN=function(k){
    m <- max.perm[max.perm$name==k & max.perm$mx < 0, c("mx", "lambda_perbase")]
    m$mx <- -1*m$mx
    m <- m[dim(m)[1]:1, ]
    get_thresh_with_rate1(m, lambda.pb, np=np)
  })
  thresh[match(keep.segs, segment.bounds$name),2] <- -1*zneg
  return(thresh)
}

#Given a per-base fdr, what threshold corresponds?
#This function assumes all peaks are positive (transform to get negative
  #thresholds)
project_thresh_with_rate1 <- function(ll, rate, np=10){
  ii <- order(c(rate-ll[,2], 0))
  N <- nrow(ll) + 1
  zero_ii <- which(ii==N)
  if((zero_ii-1) < floor(np/2)){
    nneg <- zero_ii -1
    npos <- np-nneg
  }else if((N-zero_ii) < ceiling(np/2)){
    npos <- N-zero_ii
    nneg <- np-npos
  }else{
    nneg <- floor(np/2)
    npos <- ceiling(np/2)
  }
  if(nneg == 0){
    ix <- ii[2:(np+1)]
  }else if(npos==0){
    ix <- ii[(N-np):(N-1)]
  }else{
    ix <- ii[(zero_ii -nneg):(zero_ii-1)]
    ix <- c(ix, ii[(zero_ii + 1):(zero_ii + npos)])
  }
  #ix <- order(abs(rate-ll[,2]))[1:np]
  #Doing the regression in the same way as get_rate_with_thresh makes
    #the two functions exactly reversible
  ff <- lm(log10(ll[ix, 2])~ll[ix, 1])
  slope <- 1/ff$coefficients[2]
  b <-  - (ff$coefficients[1]/ff$coefficients[2])
  return(slope*log10(rate) + b)
}

#Want to just invert get_rate_with_thresh
#Note that ll is sorted so that ll[,1] decresases and ll[,2] increases
get_thresh_with_rate1 <- function(ll, rate, np=10, tol=1e-13){
  t1 <- project_thresh_with_rate1(ll, rate, np=np)
  r1 <- get_rate_with_thresh(ll, t1, np=np)
  if(abs(log10(r1)-log10(rate)) < tol) return(t1)
  if(rate < ll[1,2]){ #Rate is smaller than any rate in ll
    d <- t1-ll[1,1]
    t <- seq(ll[2,1], t1 + (d/2), length.out=1000)
  }else if (rate < ll[2,2]){
    t <- seq(ll[3,1], ll[1,1], length.out=1000)
  }else{
    ii <- order(c(rate-ll[,2], 0))
    N <- nrow(ll) + 1
    zero_ii <- which(ii==N)
    left <- ii[max(1, zero_ii-2)]
    right <- ii[zero_ii + 2]
    t <- seq(ll[left, 1], ll[right, 1], length.out=1000)
  }
  rr <- sapply(t, FUN=function(thresh){
    fret:::get_rate_with_thresh(ll, thresh, np=10)
  })
  t2 <- approx(x=log10(rr), y=t, xout=log10(rate))$y
  r2 <- get_rate_with_thresh(ll, t2, np=np)
  if(!abs(log10(r1)-log10(rate)) < tol) cat("Warning: threshold may be inacurate\n")
  return(t2)
}

get_discoveries <- function(max1, thresholds){
  discoveries <- matrix(nrow=nrow(max1), ncol=3)
  ix <- which(thresholds$num.disc > 0)
  chru <- unique(thresholds$chrom[ix])
  max1$tpos <- thresholds$thresh.pos[match(max1$name, thresholds$name)]
  max1$tneg <- thresholds$thresh.neg[match(max1$name, thresholds$name)]
  if(any(is.na(max1$tpos)) | any(is.na(max1$tneg))) stop("Error in get_discoveries.1.\n")
  if(any(max1$tpos[max1$mx > 0] > max1$mx[max1$mx > 0]) |
            any(max1$tpos[max1$mx < 0] < max1$mx[max1$mx < 0])){
    stop("Error in get_discoveries.2.\n")
   }
  for(c in chru){
    cat(c, ": ")
    file <- unique(thresholds$file[ix][thresholds$chrom[ix]==c])
    stopifnot(length(file)==1)
    stats <- getobj(file)
    ixc <- which(max1$chr==c)
    discoveries[ixc, 3] <- c
    for(i in ixc){
      cat(i, " ")
      strt <- which(stats$sts.smooth$pos==max1$start[i])
      stp <- which(stats$sts.smooth$pos==max1$stop[i])
      D <- stats$sts.smooth[strt:stp,]
      tpos <- max1$tpos[i]
      tneg <- max1$tneg[i]
      excr <- fret:::excursions(D$ys, c(tpos, tneg))
      discoveries[i,1] <- D$pos[min(excr[,1])]
      discoveries[i,2] <- D$pos[max(excr[,2])]
    }
    cat("\n")
  }
  return(discoveries)
}
