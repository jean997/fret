

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

  if(all(m < zmin)){
    mx <- cbind(zmin, 0)
  }else if(sum(m >=zmin) < 5){
    mx <- mx[1:5, ]
  }else{
    mx <- mx[m >= zmin,]
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




#' Find thresholds for a range of lambda values. Calculate FDR.
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param file.list List of files. Each file should contain a list of objects produced by
#'discopony_maxes1
#'@param zmin Lower bound on significance thresholds.
#'@param nlam Number of lambda values to consider
#' @return A list with items z, Robs, and fdr.
#'@export
fret_choose_z <- function(file.list, zmin, nlam, log.lambda.min=NULL, log.lambda.max=NULL){
  if(is.null(log.lambda.min)){
    n.chunk=0
    log.lambda.min <- Inf
    log.lambda.max <- -Inf
    for(f in file.list){
      R <- getobj(f)
      for(i in 1:length(R)){
        log.lambda.min <- min(log.lambda.min, log10(R[[i]]$mx[,2][R[[i]]$mx[,2] > 0 &  R[[i]]$mx[,1] >=zmin]))
        log.lambda.max <- max(log.lambda.max, log10(R[[i]]$mx[,2][R[[i]]$mx[,1] >=zmin]))
      }
      n.chunk = n.chunk + length(R)
    }
    lams <- seq(log.lambda.min, log.lambda.max, length.out=nlam)
  }else{
    n.chunk=0
    for(f in file.list){
      R <- getobj(f)
      n.chunk = n.chunk + length(R)
    }
    lams <- seq(log.lambda.min, log.lambda.max, length.out=nlam)
  }

  names <- c()
  nbp <- 0
  n.seg <- length(file.list)
  z <- matrix(nrow=nlam, ncol=n.chunk + 1)
  Robs <- matrix(nrow=nlam, ncol=n.chunk+1)
  z[,1] <- Robs[, 1] <- lams
  ct <- 1
  for(f in file.list){
    R <- getobj(f)
    for(i in 1:length(R)){
      if(nrow(R[[i]]$mx)==1 & R[[i]]$mx[1, 2]==0){
        #No permutation peaks ever reached zmin --> zmin is the threshold for all lambda
        z[, ct+1] = zmin
        Robs[, ct+1] = sum(R[[i]]$max1 >= zmin)
      }else{
        z[, ct+1] <- approx(y=R[[i]]$mx[,1], x=log10(R[[i]]$mx[,2]),
                              xout=lams, yright=zmin, yleft=Inf)$y
        z[, ct+1] <- pmax(zmin, z[, ct+1])
        Robs[, ct+1] <- sapply(z[, ct+1], FUN=function(zz){sum(R[[i]]$max1 > zz)})
      }
      names <- c(names, R[[i]]$file)
      nbp <- nbp + R[[i]]$nbp
      ct = ct + 1
    }
  }
  z <- data.frame(z)
  names(z)<- c("lambda", names)
  Robs <- data.frame(Robs)
  names(Robs) <- c("lambda", names)
  fdr <- (10^(Robs[,1]))*nbp/rowSums(Robs[, -1, drop=FALSE])
  ret <- list("Robs"=Robs, "z"=z, "fdr"=fdr)
  return(ret)
}

#This function assumes a v. specific file name structure
#'@export
discopony_pull_regions <- function(results.file, thresh){
  Z <- getobj(results.file)
  if(all(Z$fdr > thresh)){
    cat("No intervals achieve significance")
    return(0)
  }
  ix <- max(which(Z$fdr <= thresh))
  R <- as.numeric(Z$Robs[ix,-1])
  z <- as.numeric(Z$z[ix, -1])
  file.names = names(Z$Robs)[-1]

  df <- data.frame(cbind(z[R > 0],  R[R > 0]), stringsAsFactors = FALSE)
  names(df) <- c("z", "R")
  df$file <- file.names[R > 0]

  ivls <- matrix(nrow=sum(df$R), ncol=5)
  #Chr chunk start stop z
  cat(sum(df$R), " intervals in ", nrow(df), " chunks", "\n")

  df$chr <- strsplit_helper(df$file, "/", 1, fixed=TRUE)
  chrnum <- as.numeric(strsplit_helper(df$chr, "chr", 2))
  df$chunk <- strsplit_helper(df$file, ".txt", 1, fixed=TRUE )
  df$chunk <- strsplit_helper(df$chunk, "chunk", 2, fixed=TRUE)
  o <- order(chrnum, as.numeric(df$chunk))
  df <- df[o, ]

  #One row per interval
  ivls[,1] <- rep(df$chr, df$R)
  ivls[,2] <- rep(df$chunk, df$R)
  ivls[,5] <- rep(df$z, df$R)
  j <- 1
  c <- ""
  for(i in 1:nrow(df)){
    if(df$chr[i]!=c){
      c <- df$chr[i]
      cat(c, "..")
      fn <- paste0("discopony_output/", c, "/", c, "_maxes.upd.RData")
      ff <- getobj(fn)
      my.names <- sapply(ff, FUN=function(f){f$file})
    }
    k <- which(my.names==df$file[i])
    iv <- name_clusters_merged(x=ff[[k]]$ys, z=df$z[i], z0=0.3*0.9)
    nn <- nrow(iv)
    stopifnot(nn==df$R[i])
    iv <- matrix(ff[[k]]$pos[as.matrix(iv)], ncol=2, byrow=FALSE)
    iv[,2] <- iv[,2]+1
    ivls[j:(j+nn-1), c(3, 4)] <- iv
    j <- j + nn
  }

  ivls <- data.frame(ivls, stringsAsFactors=FALSE)
  names(ivls) <- c("chr", "chunk", "start", "stop", "z")
  for(j in 2:5) ivls[,j] <- as.numeric(ivls[,j])
  ivls$length <- ivls$stop-ivls$start
  return(list("df"=df, "ivls"=ivls))
}

