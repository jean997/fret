#' Calculate per-base error rates for original and permutation peaks
#'@param fret.obj Either a fret object or a file containing a fret object
#'@param segment.bounds Optional data frame of segment bounds. If ommitted
#'then fret.obj should contain an item `seg.bounds` which will be used instead.
#'@param parallel Run in parallel?
#'@export
fret_rates_prelim <- function(fret.obj, segment.bounds=NULL,
                              parallel=FALSE, save.file=NULL){
  if(class(fret.obj)=="character") fret.obj <- getobj(fret.obj)

  stopifnot(all(c("m1", "mperm") %in% names(fret.obj)))

  if(is.null(segment.bounds)){
    stopifnot("seg.bounds" %in% names(fret.obj))
    segment.bounds <- fret.obj$seg.bounds
  }
  stopifnot(ncol(segment.bounds)==3)
  stopifnot(names(segment.bounds)==c("chr", "start", "stop"))
  s <- length(fret.obj$zmin)
  stopifnot(s %in% c(1, 2))
  K <- nrow(segment.bounds)
  cat("There are ", K, " segments.\n")
  segment.bounds$nbp <- segment.bounds$stop-segment.bounds$start + 1
  if(s==1){
    segment.bounds$max_lambda_perbase <- rep(NA, K)
  }else{
    segment.bounds$max_lambda_perbase_pos <- rep(NA, K)
    segment.bounds$max_lambda_perbase_neg <- rep(NA, K)
  }
  mlp_ix <- grep("max_lambda", names(segment.bounds))

  #Match peaks to segments

  fret.obj$m1$segment <- match_segments(chr=fret.obj$m1$chr, pos=fret.obj$m1$pos, segment.bounds = segment.bounds, parallel=parallel)

  nsegs1 <- sapply(1:K, FUN=function(i){ sum(fret.obj$m1$segment==i)})

  fret.obj$mperm$segment <- match_segments(chr=fret.obj$mperm$chr, pos=fret.obj$mperm$pos, segment.bounds = segment.bounds, parallel=parallel)

  nsegsperm <- sapply(1:K, FUN=function(i){ sum(fret.obj$mperm$segment==i)})

  #For each peak in m1 and mperm we want to find the value of lambda
  #at the max height of the peak.
  #We also want to find the overall maximum achievable lambda in each segement.
  #This is important since some segments have a very low max lambda
  #(e.g. we almost never see a permutation peak above z0 in these segments)

  fret.obj$m1$lambda_perbase <- rep(0, nrow(fret.obj$m1))
  fret.obj$mperm$lambda_perbase <- rep(0, nrow(fret.obj$mperm))
  for(i in 1:K){
    if(i %% 10 == 1) cat(i, "..")
    #m1.ix <- which(max1$segment==i)
    #perm.ix <- which(max.perm$segment==i)
    #No peaks in data and almost no peaks in permutations
    if(nsegs1[i]==0 & nsegsperm[i] < 2){
      segment.bounds[i, mlp_ix] <- 0
      next
    }
    #No permutation peaks above z0 but there are peaks above zmin in data
      #i.e. highly significant but cant estimate significance
      #because all perm peaks are too low
      #probably very rare
    if(nsegs1[i] > 0 & nsegsperm[i] ==0){
      m1.ix <- which(max1$segment==i)
      segment.bounds[i, mlp_ix] <- 0
      fret.obj$m1$lambda_perbase[m1.ix] <- 0
      next
    }
    #Rates for permutation peaks
    perm.ix <- which(fret.obj$mperm$segment==i)
    m <- fret.obj$mperm$mx[perm.ix]
    o <- order(m, decreasing=TRUE)
    oinv <- match(1:length(m), o)
    ll <- fret:::lamtab(mx=m, zmin=fret.obj$zmin, nbp = segment.bounds$nbp[i], n.perm=fret.obj$n.perm)
    fret.obj$mperm$lambda_perbase[perm.ix] <- ll[,2][oinv]
    if(s==1){
      segment.bounds[i, mlp_ix] <- fret:::get_rate_with_thresh(ll, fret.obj$zmin, np=2)
    }else{
      segment.bounds[i, mlp_ix[1]] <- fret:::get_rate_with_thresh(ll, fret.obj$zmin[1], np=2)
      segment.bounds[i, mlp_ix[2]] <- fret:::get_rate_with_thresh(ll, fret.obj$zmin[2], np=2)
    }
    fret.obj$mperm[perm.ix, ] <- fret.obj$mperm[perm.ix,][o,]
    if(nsegs1[i] > 0){
      m1.ix <- which(fret.obj$m1$segment==i)
      rts <- sapply(fret.obj$m1$mx[m1.ix], FUN=function(thresh){
          fret:::get_rate_with_thresh(ll, thresh)
      })
      fret.obj$m1$lambda_perbase[m1.ix] <- rts
    }
  }
  fret.obj$seg.bounds <- segment.bounds
  if(!is.null(save.file)){
    save(fret.obj, file=save.file)
    return(0)
  }
  return(fret.obj)
}


#In one segment, at a partiucular threshold, what is the rate of false discoveries?
#ll produced by lamtab. Two columns. First column is threshold. Second column is perbase fdr
#ll is assumed to be sorted so that the first column is decreasing.
#np If thresh is larger than the largest threshold or smaller than the smallest threshold
# in the first column of ll then we estimate the rate by fitting a line and projecting.
#np controls how many points to use when fitting the line.
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