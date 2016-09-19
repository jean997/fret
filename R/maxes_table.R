#'@export
maxes_table <- function(dat.file, pheno.file, s0, zmin, seed, n.perm,
                        z0=zmin*0.3, save.perm.stats=FALSE, range=NULL,
                        bandwidth=50, out.file=NULL, chrom="chr1",
                        stat.func=huber_stats, ...){

  stopifnot(length(z0)==length(zmin))
  s <- length(zmin)
  stopifnot(s %in% c(1, 2))
  name.root <- unlist(strsplit(dat.file, ".txt"))[1]
  if(is.null(out.file)) out.file <- paste0(name.root, "_out.RData")

  #Read phenotype
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
  dat <- read_delim(dat.file, delim=" ")
  if(nrow(dat)==0){
    cat("No data in ", dat.file)
    return(0)
  }
  if(! names(dat)[1]=="pos") stop("First column of data matrix should be genomic position and be named 'pos'")
  pos <- dat$pos
  #Make sure the phenotype is sorted correctly
  if(!all(names(dat)[-1] %in% X[[1]])) stop("Not all of the samples in the data file are in the phenotype file.\n")
  X <- X[match(names(dat)[-1], X[[1]]),  ]

  #Range
  if(!is.null(range)){
    stopifnot(length(range)==2)

    ixstart <- max(1, min(which(pos >=(range[1]-ceiling(bandwidth/2)) )))
    ixstop <- min(length(pos), max(which(pos <=(range[2]+ceiling(bandwidth/2)) )))
    dat <- dat[ixstart:ixstop,]
    pos <- pos[ixstart:ixstop]

    ix1 <- min(which(pos >=range[1]))
    ix2 <- max(which(pos <= range[2]))
  }else{
    ix1 <- 1
    ix2 <- length(pos)
    range <- c(pos[1], pos[length(pos)])
  }

  pos.out = pos[ix1:ix2]


  #Calculate statistics
  cat("Calculating test statistics..\n")
  sts <- stat.func(Y=dat[, -1], x=X[[2]],s0=s0, ...)
  ys <- ksmooth_0(x=pos, y=sts[3,], bandwidth = bandwidth)[ix1:ix2]
  sts <- data.frame(t(sts))
  names(sts) <- c("Beta", "SD", "stat")
  sts$pos <- pos
  #Smooth statistics
  stats.smooth <- data.frame("pos"=pos.out, "ys"=ys)


  #Find peak heights
  #Intervals defined by z0:
  max1 <- mxlist(ys, z0, zmin, return.ix = TRUE)
  if(s==1) max1 <- max1[max1$mx > zmin,]
    else max1 <- max1[max1$mx > zmin[1] | max1$mx < zmin[2],]
  if(nrow(max1) == 0){
    cat("No clusters exceed ", zmin, "\n")
  }
  #Convert index to position
  max1$pos <- pos.out[max1$ix]
  max1$iv1 <- pos.out[max1$iv1]
  max1$iv2 <- pos.out[max1$iv2]
  max1$chr <- chrom
  max1 <- max1[, c("mx", "chr", "pos", "iv1", "iv2")]
  if(n.perm==0){
    R <- list("max1"=max1,"file"=dat.file,
              "stats" = sts, "stats.smooth"=stats.smooth,
              "z0"=z0, "zmin"=zmin)
    return(R)
  }

  cat("Calculating permutation statistics...\n")

  max.perm <- data.frame("mx"=c(), "ix"=c(), "iv1"=c(), "iv2"=c())
  stats.perm <- apply(perms, MARGIN=2, FUN=function(l){
    huber_stats(dat[,-1], x=l, s0=s0, maxit=maxit)[3,]
  })
  stats.perm.smooth <- apply(stats.perm, MARGIN=2, FUN=function(yy){
    ksmooth_0(x=pos, y=yy, bandwidth = bandwidth)[ix1:ix2]
  })
  vv <- apply(stats.perm.smooth, MARGIN=1, FUN=var)

  max.perm <- apply(stats.perm.smooth, MARGIN=2, FUN=function(yys){
    mxlist(yys, z0, zmin, return.ix = TRUE)
  })
  max.perm <- do.call(rbind, max.perm)

  if(nrow(max.perm) > 0){
    #Convert positions
    max.perm$pos <- pos.out[max.perm$ix]
    max.perm$iv1 <- pos.out[max.perm$iv1]
    max.perm$iv2 <- pos.out[max.perm$iv2]
    max.perm$chr <- chrom
    max.perm <- max.perm[order(max.perm$pos, decreasing = FALSE), c("mx", "chr", "pos", "iv1", "iv2")]
  }
  if(save.perm.stats){
    R <- list("max1"=max1, "max.perm"=max.perm, "file"=dat.file,
              "z0"=z0, "zmin"=zmin, "n.perm"=n.perm, "perm.var"=vv,
              "stats.smooth"=stats.smooth,"stats.perm" = stats.perm,
              "stats.perm.smooth"= stats.perm.smooth)
  }else{
    R <- list("max1"=max1, "max.perm"=max.perm, "file"=dat.file,
              "z0"=z0, "zmin"=zmin, "n.perm"=n.perm, "perm.var"=vv,
              "stats.smooth" = stats.smooth)
  }
  save(R, file=out.file)
  return(nrow(max.perm))
}


#Don't save everything
#if(s==1){
#  ix <- which(max.perm$mx > zmin)
#  if(length(ix) < 10){
#    n <- min(10, nrow(max.perm))
#    ix <- order(max.perm$mx, decreasing=TRUE)[1:n]
#  }
#}else{
#  ixpos <- which(max.perm$mx > zmin[1])
#  if(length(ixpos) < 10){
#    n <- min(10, sum(max.perm$mx > 0))
#    ixpos <- order(max.perm$mx, decreasing=TRUE)[1:n]
#  }
#  ixneg <- which(max.perm$mx < zmin[2])
#  if(length(ixneg) < 10){
#    n <- min(10, sum(max.perm$mx < 0))
#    ixneg <- order(max.perm$mx, decreasing=FALSE)[1:n]
#  }
#  ix <- c(ixpos, ixneg)
#}
#max.perm <- max.perm[ix,]
