

#' @title Calculate FRET test statistics
#' @param pheno.file File containing genomic phenotype
#' @param trait.file File containing trait information
#' @param s0 Variance inflation constant
#' @param seed Seed
#' @param n.perm Number of permutations.
#' If n.perm=0, only test statistics for unpermuted data will be calculated
#' @param zmin Minimum threshold (may be missing if n.perm=0)
#' @param z0 Merging threshold. If missing z0 = 0.3*zmin.
#' @param pheno.transformation If the phenotype is to be transformed, provide a function taking one argument
#' @param trait Name of trait (should match header of trait.file)
#' @param covariates List of covariates found in trait.file to adjust for
#' @param save.perm.stats Should the permutation test statistics be saved?
#' @param range Base-pair range. If missing, stats for the whole file will be calculated
#' @param chunksize Size of chunks to read at one time
#' @param bandwidth Smoothing bandwidth
#' @param smooother Either "ksmooth_0" or "ksmooth". The ksmooth_0 smoother treats base-pairs
#' with no data as though the phenotype at that position is equal to zero for all samples (appropriate for DNase-seq data).
#' The ksmooth smoother treats these base-pairs as missing (appropriate for bisulfite sequencing data).
#' @param stat.type Either "huber" or "lm"
#' @param maxit Maximum iterations for Huber estimator.
#' @param out.file Name an output file
#' @param chrom Name of chromosome to print in output file
#'@export
fret_stats3 <- function(pheno.file, trait.file, s0, seed, n.perm, zmin=NULL, z0=zmin*0.3,
                        pheno.transformation=NULL, trait=c("x"), covariates=c(),
                        range=NULL, chunksize=1e5,
                        bandwidth=50, smoother=c("ksmooth_0", "ksmooth", "none"),
                        stat.type=c("huber", "lm"), maxit=50,
                       out.file=NULL, chrom="chr1", parallel=TRUE, tmp.dir ="/tmp/"){
  ############
  #  Options #
  ############
  if(length(trait) !=1) stop("ERROR: Handling of multivariate traits is not implemented yet!\n")
  if(is.null(zmin)) stopifnot(n.perm==0)
  if(smoother=="none") stopifnot(n.perm==0)
  #stat type
  stat.type <- match.arg(stat.type)
  if(stat.type=="huber"){
    lm.func <- function(formula, data){
      rlm(formula, data=data, maxit=maxit)
    }
    if(!parallel) stat.func <- function(Y, x, s0){ huber_stats(Y, x, s0, maxit=maxit)}
    	else stat.func <- function(Y, x, s0){ huber_stats_parallel(Y, x, s0=s0, maxit=maxit)}
  }else if(staty.type=="lm"){
    lm.func <- function(formula, data){
      lm(formula, data=data)
    }
    stat.func <- function(Y, x, s0){ lm_stats_parallel(Y, x, s0=s0)}
  }
  #smoother type
  smoother <- match.arg(smoother)
  if(smoother=="ksmooth_0"){
    smooth.func <- function(x, y, xout, bandwidth){ksmooth_0(x, y,xout, bandwidth, stitch=1e5, parallel=parallel)}
  }else if(smoother=="ksmooth"){
    smooth.func <- function(x, y, bandwidth){ ksmooth(x=x, y=y, x.points=x, bandwidth=bandwidth)$y}
  }

  stopifnot(length(z0)==length(zmin) | is.null(zmin))
  s <- length(zmin)
  stopifnot(s %in% c(0, 1, 2))
  if(!is.null(pheno.transformation)) stopifnot("function" %in% class(pheno.transformation))
  chunksize <- floor(chunksize)

  #Object to return
  R <- list("pheno.file"=pheno.file, "trait.file"= trait.file,
            "range"=range, "trait"=trait, "covariates"=covariates,
            "pheno.transformation"=pheno.transformation,
            "bandwidth"=bandwidth, "smoother"=smoother,
            "z0"=z0, "zmin"=zmin)


  ####################
  #  Read trait data #
  ####################
  X <- read_delim(trait.file, delim=" ")
  if(!all(trait %in% names(X))) stop("ERROR: I didn't find colunns matching the specified trait name in the phenotype file.\n")
  if(!all(covariates %in% names(X))) stop("ERROR: I didn't find colunns matching the specified covariates in the phenotype file.\n")
  n <- nrow(X)
  #Adjust trait for covariates
  if(length(covariates) > 0){
    ff <- as.formula(paste0(trait, "~", paste0(covariates, collapse="+")))
    fitx <- lm.func(ff, X)
    x <- fitx$residuals
  }else{
    x <- X[[trait]]
  }
  #Permuted phenotypes
  if(n.perm > 0){
    set.seed(seed)
    perms <- replicate(n=n.perm, expr = {
      sample(x, size=n, replace=FALSE)
    })
  }

  ###################
  # Read pheno data #
  ###################
  h <- readLines(pheno.file, n=1)
  h <- unlist(strsplit(h, " "))

  #Make sure the trait data is sorted correctly
  if(!all(h[-1] %in% X[[1]])) stop("Not all of the samples in ", pheno.file, " are in ", trait.file, ".\n")
  x <- x[match(h[-1], X[[1]])]
  X <- X[match(h[-1], X[[1]]), ]

  #Use LaF package to read file in chunks
  dm <- detect_dm_csv(filename=pheno.file, header=TRUE, sep=" ")
  df.laf <- laf_open(dm)

  nl <- as.numeric(strsplit(system(paste0("wc -l ", pheno.file), intern=TRUE), " ")[[1]][1])
  nchunk <- ceiling(nl/chunksize)
  cat(pheno.file, " contains ", nl, " lines. I will process it in ", nchunk, " chunks.\n")
  #Adjust Range
  if(!is.null(range)){
    #We need to look at data beyond the specified range to get the smoothing right
    #And also because some peaks might go outside of the range
    new.range <- c(range[1]-3*bandwidth , range[2]+3*bandwidth)
    new.range[1] <- max(1, new.range[1])
  }

  done <- FALSE
  read.ct <- 1
  keep.ct <- 1
  #For saving as we go so we don't use an absurd amount of memory
  tdir <- paste0(tmp.dir, paste0(sample(c(letters, LETTERS), size=4), collapse=""), "_", chrom)
  system(paste0("mkdir ", tdir))

  while(!done){
    cat("Chunk ", read.ct, "(", keep.ct, ")\n")
    #read data
    goto(df.laf, max(0, (read.ct-1)*chunksize-bandwidth))
    dat <- next_block(df.laf, nrows=chunksize + 2*bandwidth)
    cat("dat: ", dim(dat), "\n")
    if(nrow(dat) < chunksize + 2*bandwidth) done <- TRUE
    if(nrow(dat) <= bandwidth) break

    if(read.ct==1) nstart <- 1
      else nstart <- bandwidth + 1
    nend <- min(nrow(dat), nstart + chunksize-1)
    if(!is.null(range)){
      if(dat[[1]][nend] < new.range[1]){
        read.ct <- read.ct + 1
        next
      }
      if(dat[[1]][nend] > new.range[2]){
        done <- TRUE
        nend <- which.min(dat[[1]] < new.range[2])
      }
      if(dat[[1]][nstart] < new.range[1] & dat[[1]][nend] > new.range[1]){
        nstart <- which.max(dat[[1]] >= new.range[1])
      }
    }

    ###################################
    # Adjust phenotype for covariates #
    ###################################
    if(!is.null(pheno.transformation)){
      cat("Adjusting phenotype.\n")
      dat[,-1] <- apply(dat[,-1], MARGIN=2, FUN=function(y){
        pheno.transformation(y)
      })
      if(any(is.na(dat[,-1])) | any(!is.finite(dat[,-1]))) stop("ERROR: Phenotype adjustment gives infinite or missing value.\n")
    }
    if(length(covariates) > 0){
      cat("Regressing covariates on phenotype.\n")
      dat[,-1] <- apply(dat[,-1], MARGIN=1, FUN=function(y){
        ff <- as.formula(paste0("y~", paste0(covariates, collapse="+")))
        fity <- lm.func(ff, X)
        fity$residuals
      })
    }

    ###############################################
    # Calculate (non-permutation) test statistics #
    ###############################################
    cat("Calculating test statistics..\n")
    my.sts <- stat.func(Y=dat[, -1], x=x,s0=s0)
    my.sts <- data.frame(t(my.sts))
    names(my.sts) <- c("Beta", "SD", "stat")
    my.sts$pos <- dat[[1]]
    if(keep.ct ==1 ) sts <- my.sts[nstart:nend,]
      else sts <- rbind(sts, my.sts[nstart:nend,])

    if(smoother=="none"){
      read.ct <- read.ct + 1
      keep.ct <- keep.ct + 1
      next
    }

    ######################################
    # Smooth (non-permutation) statistics#
    ######################################
    cat("Smoothing..\n")
    ys <- smooth.func(x=dat[[1]], y=sts$stat,xout=dat[[1]][nstart:nend], bandwidth = bandwidth)
    if(keep.ct == 1) sts.smooth <- data.frame("pos"=dat[[1]], "ys"=ys[ix1:ix2])
      else sts.smooth <- cbind(sts.smooth, data.frame("pos"=dat[[1]], "ys"=ys[ix1:ix2]))

    if(n.perm==0){
      read.ct <- read.ct + 1
      keep.ct <- keep.ct + 1
      next
    }

    #########################################
    # Calculate permutation test statistics #
    #########################################
    cat("Calculating permutation statistics...\n")
    sts.perm <- apply(perms, MARGIN=2, FUN=function(l){
      stat.func(dat[,-1], x=l, s0=s0)
    })
    ######################################
    # Smooth permutation test statistics #
    ######################################
    cat("Smoothing permutation statistics...\n")
    stat.ix <- seq(3, nrow(sts.perm), by=3)
    stats.perm.smooth <- apply(sts.perm[ix,], MARGIN=2, FUN=function(yy){
      smooth.func(x=pos, y=yy, bandwidth = bandwidth)
    })
    save(sts.perm.smooth, file=paste0(tdir, "/smooth_stats_", keep.ct, ".RData"))

    read.ct <- read.ct + 1
    keep.ct <- keep.ct + 1
  }
  rm(dat)

  R$stats <- sts
  if(smoother=="none"){
    if(!is.null(out.file)){
      save(R, file=out.file)
      return(NULL)
    }
    return(R)
  }

  R$stats.smooth <- sts.smooth
  if(is.null(zmin)){
    if(!is.null(out.file)){
      save(R, file=out.file)
      return(NULL)
    }
    return(R)
  }

  #######################################
  # Find (non-permutation) peak heights #
  #######################################
  max1 <- mxlist(sts.smooth$ys, z0, zmin, return.ix = TRUE)
  if(s==1) max1 <- max1[max1$mx > zmin,] #Unsigned
    else max1 <- max1[max1$mx > zmin[1] | max1$mx < zmin[2],] #signed
  #Convert index to position
  max1$pos <- pos[max1$ix]
  max1$iv1 <- pos[max1$iv1]
  max1$iv2 <- pos[max1$iv2]
  max1$chr <- chrom
  max1 <- max1[, c("mx", "chr", "pos", "iv1", "iv2")]
  if(!is.null(range)) max1 <- max1[max1$mx >= range[1] & max1$mx <= range[2],]

  R$max1 <- max1
  if(n.perm==0){
    if(!is.null(out.file)){
      save(R, file=out.file)
      return(NULL)
    }
    return(R)
  }

  cat("Collecting permutation test stats.\n")
  sts.perm.smooth <- matrix(nrow=0, ncol=2)
  for(i in 1:keep.ct-1){
    ss <- getobj(paste0(tdir, "/smooth_stats_", i, ".RData"))
    sts.perm.smooth <- rbind(sts.perm.smooth, ss)
    unlink(paste0(tdir, "/smooth_stats_", i, ".RData"))
  }
  R$perm.var <- apply(stats.perm.smooth, MARGIN=1, FUN=var)
  max.perm <- apply(stats.perm.smooth, MARGIN=2, FUN=function(yys){
    mxlist(yys, z0, zmin, return.ix = TRUE)
  })
  R$max.perm <- do.call(rbind, max.perm)


  if(!is.null(out.file)){
    save(R, file=out.file)
    return(nrow(max.perm))
  }
  return(R)
}
