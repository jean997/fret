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

  tot.disc <- sum(obj$max1$fdr <= target.fdr) ###This is the number of discoveries
  #We want to draw thresholds with lambda = target.fdr*total num discoveries
  lam.target <- target.fdr*tot.disc
  stopifnot(sum(obj$max1$lambda <= lam.target)==sum(obj$max$fdr <= target.fdr))
  segs <- unique(obj$max1$name[obj$max1$fdr <= target.fdr])
  tt <- fret:::get_thresh_with_rate(obj$max.perm, obj$segment.bounds,
                                    lam.target, obj$zmin, np=2, segs=segs)
  thresholds$thresh.pos <- tt[,1]
  if(s==1) thresholds$thresh.neg <- -tt[,1]
  else thresholds$thresh.neg <- tt[,2]

  for(j in 1:K) thresholds$num.disc[j] <- sum(obj$max1$name[1:tot.disc]==obj$segment.bounds$name[j])
  thresholds$name <- obj$segment.bounds$name
  ix <- which(thresholds$num.disc > 0)
  thresholds$chrom[ix] <- obj$max1$chr[match(ix, obj$max1$segment)]
  if(is.null(stats.files)) return(thresholds)
  cat("Thresholds found. Retrieving discoveries.\n")
  thresholds$file[ix] <- stats.files[match(thresholds$chrom[ix], stats.files[,1]), 2]
  discoveries <- get_discoveries(max1=obj$max1[1:tot.disc,], thresholds = thresholds)
  ret <- list("thresholds"=thresholds, "discoveries"=discoveries)
  return(ret)
}

get_thresh_with_rate <- function(max.perm, segment.bounds,
                                 lambda, zmin, np=2, segs=NULL){
  s <- length(zmin)
  K <- nrow(segment.bounds)
  stopifnot("name" %in% names(segment.bounds))
  stopifnot("nbp" %in% names(segment.bounds))
  mlp_ix <- grep("max_lambda_perbase", names(segment.bounds))
  stopifnot(length(mlp_ix) > 0)

  #stopifnot(all(dim(max.lambda.pb)==dim(nbp)))
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
get_thresh_with_rate1 <- function(ll, rate, np=4){
  if(rate < min(ll[,2])){
    ff <- lm(ll[1:np, 1]~log10(ll[1:np, 2]))
    return(ff$coefficients[2]*log10(rate) + ff$coefficients[1])
  }
  return(approx(y=ll[,1], x=log10(ll[,2]),
                xout=log10(rate), yright=min(ll[,1]))$y)
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
