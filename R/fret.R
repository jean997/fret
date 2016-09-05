

#' Get peak heights and permuted peak heights for one segment
#'@description Finds peak heights and permuatation peak heights for one segment
#'@param dat.file Data file. Should have first
#'column as position and a header giving sample names. One row per base pair.
#'@param pheno.file Phenotype file. First column should correspond to the header of dat.file.
#'@param s0 Variance inflation constant.
#'@param zmin Lower bound on significance thresholds.
#'@param seed Set a seed for running permutations (Required)
#'@param n.perm Number of permutations
#'@param save.perm.stats Do you want to save all the permutation statistics?
#'@param segment.bounds Length 2 vector giving start and stop base-pair for the segment.
#'If NULL, all base-pairs in dat.file are used.
#'@param out.file Output file (optional). If NULL will append "_mx.RData" to the leading portion of the
#'dat.file name.
#'@param max1.file File with data from running with n.perm=0
#'@param z0 Reference level for merging
#'@param bandwidth Smoother bandwidth.
#'@param maxit Maximum iterations for rlm.
#' @return Saves an RData file with maxes list and other info. Returns the length of that list.
#'@export
maxes_1seg <- function(dat.file, pheno.file, s0, zmin,
                      seed, n.perm, save.perm.stats=FALSE,
                      segment.bounds = NULL,
                      out.file=NULL, max1.file=NULL,
                      z0=zmin*0.3, bandwidth=50, maxit=50){

  stopifnot(length(z0)==length(zmin))
  stopifnot(length(z0) %in% c(1, 2))


  X <- read_delim(pheno.file, col_names=FALSE, delim=" ")
  n <- nrow(X)

  #Permuted phenotypes
  if(n.perm > 0){
    set.seed(seed)
    perms <- replicate(n=n.perm, expr = {
      sample( X[[2]], size=n, replace=FALSE)
    })
  }

  #Read data
  name.root <- unlist(strsplit(dat.file, ".txt"))[1]
  dat <- read_delim(dat.file, delim=" ")
  if(nrow(dat)==0) return(0)
  if(! names(dat)[1]=="pos") stop("First column of data matrix should be named 'pos'")
  pos <- dat$pos
  #Make sure the phenotype is sorted correctly
  if(!all(names(dat)[-1] %in% X[[1]])) stop("Not all of the samples in the data file are in the phenotype file.\n")
  X <- X[match(names(dat)[-1], X[[1]]),  ]

  #Chunks will have a little extra data to get the smoothing right
  if(!is.null(segment.bounds)){
    stopifnot(length(segment.bounds)==2)
    ix1 <- min(which(pos >=segment.bounds[1]))
    ix2 <- max(which(pos <= segment.bounds[2]))
  }else{
    ix1 <- 1
    ix2 <- length(pos)
    segment.bounds <- c(pos[1], pos[length(pos)])
  }
  len <- segment.bounds[2]-segment.bounds[1]+ 1
  pos.out = pos[ix1:ix2]

  #Calculate statistics
  if(is.null(max1.file)){
    sts <- huber_stats(Y=dat[, -1], x=X[[2]],s0=s0, maxit=maxit)
    ys <- ksmooth_0(x=pos, y=sts[3,], bandwidth = bandwidth)[ix1:ix2]
    sts <- data.frame(t(sts))
    names(sts) <- c("Beta", "SD", "t-stat")
    sts$pos <- pos
  }else{
    m1 <- getobj(max1.file)
    ys <- m1$stats.smooth$ys
  }

  stats.smooth <- data.frame("pos"=pos.out, "ys"=ys)
  #Find peak heights
  #Intervals defined by z0:
  max1 <- mxlist(ys, z0, zmin)
  if(length(max1) == 0){
    cat("No clusters exceed ", zmin, "\n")
    #unlink(file.name)
    return(0)
  }

  if(n.perm==0){
    R <- list("max1"=max1,"file"=dat.file, "nbp"=len,
              "stats" = sts, "stats.smooth"=stats.smooth,
              "z0"=z0, "zmin"=zmin)
    return(R)
  }

  cat("Calculating permutation statistics and peak heights.\n")
  if(save.stats){
    max.perm <- c()
    stats.perm <- array(dim=c(length(pos), 3, n.perm))
    stats.perm.smooth <- array(dim=c(length(pos.out), n.perm))
    for(i in 1:n.perm){
      cat(i, " ")
      stats.perm[ , , i] <- t(huber_stats(dat[,-1], x=perms[,i], s0=s0, maxit=maxit))
      stats.perm.smooth[,i] <- ksmooth_0(x=pos, y=stats.perm[, 3, i], bandwidth = bandwidth)[ix1:ix2]
      max.perm <- c(max.perm, mxlist(stats.perm.smooth[,i], z0, zmin))
    }
  }else{
    max.perm <- apply(perms, MARGIN=2, FUN=function(l){
      yy <- huber_stats2(dat[,-1], labs=l, s0=s0, maxit=maxit)
      yys <- ksmooth_0(x=pos, y=yy, bandwidth = bandwidth)[ix1:ix2]
      mxlist(yys, z0, zmin)
    })
  }
  m <- sort(unlist(max.perm), decreasing=TRUE)

  mx <- cbind(m, (1:length(m))/(n.perm*len))

  if( length(m)==0){
    mx <- cbind(zmin, 0)
  }

  if(is.null(out.file)){
    out.file <- paste0(name.root, "_mx.RData")
  }
  if(save.stats){
    R <- list("max1"=max1, "mx"=mx, "file"=dat.file,
              "z0"=z0, "zmin"=zmin,
              "stats.smooth"=stats.smooth,"stats.perm" = stats.perm,
              "stats.perm.smooth"= stats.perm.smooth,
              "nbp"=len)
  }else{
    R <- list("max1"=max1, "mx"=mx, "file"=dat.file,
            "z0"=z0, "zmin"=zmin,
            "nbp"=len, "stats.smooth" = stats.smooth)
  }
  save(R, file=out.file)
  return(nrow(mx))
}
