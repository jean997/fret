#' Select threshold value in each segment for a range of lambda values
#'@description Selects a series of threshold values so that each segment has approximately the same false cluster rate
#'@param perm.stats p x B matrix of permutation statistics
#'@param nlam Number of lambda values to consider
#'@param bw Binwidth for smoothing
#'@param pos Positions
#'@param z0 Reference level for merging
#'@param max.log.lam Maximum value of log(lambda) to consider
#'@param seg.ends List of segment ends. If NULL assume there is only one segment
#'@param return.z if FALSE only return maxes list
#'@param keep.lists If return.z=TRUE also return maxes list
#'@param smooth.func Either ksmooth or ksmooth_0
#' @return If return.z =TRUE, return a matrix that is nlam x (nseg + 1). The first column is log.lambda.
#' Subsequent columns give the threshold for each segment. If keep.lists=TRUE or return.z=FALSE (also) return
#' a list of cluster peaks for each segment. mx has length nseg and each element has two columns - peak height and
#' empirical lambda value.
#'@export
choose_z_even <- function(perm.stats, nlam, bw, pos, z0, x.range=NULL,
                          max.log.lam=NULL, seg.ends=NULL, return.z=TRUE,
                          keep.lists=FALSE, smooth.func = c("ksmooth", "ksmooth_0")){
  #How many segments are there?
  if(is.null(seg.ends)){
    seg.ends <- c(nrow(perm.stats))
  }
  nseg <- length(seg.ends)
  nperm <- ncol(perm.stats)
  if(is.null(x.range)) pos.out <- pos
    else pos.out <- pos[x.range[1]:x.range[2]]

  #Smoothing function
  smooth.func <- match.arg(smooth.func)
  if(smooth.func =="ksmooth"){
    smooth.func <- function(x){
      ksmooth(x=pos, y=x, bandwidth=bw, x.points=pos)$y
    }
  }else if(smooth.func == "ksmooth_0"){
    smooth.func <- function(x){
      ksmooth_0(pos, x, bandwidth=bw)
    }
  }
  #Smooth the statistics
  perm.smooth <- apply(perm.stats, MARGIN=2, FUN=function(x){
    smooth.func(x)
  })

  strt <- 1
  mx <- list()
  for(i in 1:nseg){
    cat(i, " ")
    maxes <- apply(perm.smooth[strt:seg.ends[i],], MARGIN=2, FUN=function(xs){
      if(all(abs(xs) <= z0)){
        c()
      }else{
        q0 <-rle( abs(xs) > z0 )
        p0 <- length(q0$lengths)
        starts0 <- c(1, cumsum(q0$lengths)[-p0]+1)[q0$values]
        stops0 <- (cumsum(q0$lengths))[q0$values]
        sapply(1:length(starts0), FUN=function(j){ max(abs(xs)[starts0[j]:stops0[j]])})
      }
    })
    m <- sort(unlist(maxes), decreasing=TRUE)
    mx[[i]] <- cbind(m, (1:length(m))/(nperm*(seg.ends[i]-strt + 1)))
    strt <- seg.ends[i] + 1
  }
  cat("\n")
  if(!return.z) return(mx)


  min.log.lam <- min(sapply(mx, FUN=function(m){min(log10(m[,2]))}))

  if(is.null(max.log.lam)) max.log.lam <- min(sapply(mx, FUN=function(m){max(log10(m[,2]))}))

  lams <- seq(min.log.lam, max.log.lam, length.out=nlam)
  z <- sapply(mx, FUN=function(m){
   approx(y=m[,1], x=log10(m[,2]), xout=lams)$y
  })
  z <- data.frame(cbind(lams, z))
  names(z) <- c("lambda", paste0("z", 1:nseg))
  if(! keep.lists) return(z)

  return(list("mx"=mx, "z"=z))

}
