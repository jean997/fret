

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
#' @param chrom Name of chromosome
#'@export
fret_stats <- function(pheno.file, trait.file, s0, seed, n.perm, zmin=NULL, z0=zmin*0.3,
                        pheno.transformation=NULL, trait=c("x"), covariates=c(),
                        bandwidth=150, smoother=c("ksmooth_0", "ksmooth", "none"),
                        stat.type=c("huber", "lm"), maxit=50, margin=3*bandwidth,
                        range=NULL, chunksize=1e5, which.chunks=NULL,
                       temp.prefix=NULL, temp.dir ="./",
                       chrom="chr1", parallel=TRUE, cores=parallel::detectCores()-1){
  ############
  #  Options #
  ############
  if(length(trait) !=1) stop("ERROR: Handling of multivariate traits is not implemented yet!\n")
  if(is.null(zmin)) stopifnot(n.perm==0)

  #stat type
  stat.type <- match.arg(stat.type)
  if(stat.type=="huber"){
    lm.func <- function(formula, data){
      rlm(formula, data=data, maxit=maxit)
    }
    if(!parallel) stat.func <- function(Y, x, s0){ huber_stats(Y, x, s0, maxit=maxit)}
    	else stat.func <- function(Y, x, s0){ huber_stats_parallel(Y, x, s0=s0, maxit=maxit, cores=cores)}
  }else if(stat.type=="lm"){
    lm.func <- function(formula, data){
      lm(formula, data=data)
    }
    stat.func <- function(Y, x, s0){ lm_stats_parallel(Y, x, s0=s0, cores=cores)}
  }
  #smoother type
  smoother <- match.arg(smoother)
  if(smoother=="none") stopifnot(n.perm==0)
  if(smoother=="ksmooth_0"){
    smooth.func <- function(x, y, xout, bandwidth){
      ksmooth_0(x, y,xout, bandwidth, stitch=floor(chunksize/100), parallel=parallel, cores=cores)
    }
  }else if(smoother=="ksmooth"){
    smooth.func <- function(x, y, xout, bandwidth){
      ksmooth(x=x, y=y, x.points=xout, bandwidth=bandwidth)$y
    }
  }
  #Signed or unsigned?
  stopifnot(length(z0)==length(zmin) | is.null(zmin))
  s <- length(zmin)
  stopifnot(s %in% c(0, 1, 2))
  if(!is.null(pheno.transformation)) stopifnot("function" %in% class(pheno.transformation))
  #Chunks
  chunksize <- floor(chunksize)
  #nl <- as.numeric(strsplit(system(paste0("wc -l ", pheno.file), intern=TRUE), " ")[[1]][1]) -1
  nl <- determine_nlines(pheno.file)-1
  nchunk <- max(1, ceiling(nl/chunksize))
  cat(pheno.file, " contains ", nl, " lines. It will be processed it in ", nchunk, " chunks.\n")
  if(!is.null(which.chunks)){
    stopifnot(chunksize > 0)
    if(any(which.chunks > nchunk)) stop("ERROR: Some requested chunks are too large")
    cat("This job will analyze ", length(which.chunks), " chunks.\n")
  }else{
    which.chunks <- 1:nchunk
  }
  #Adjust Range
  if(!is.null(range)){
    #We need to look at data beyond the specified range to get the smoothing right
    #And also because some peaks might go outside of the range
    new.range <- c(range[1]-margin , range[2]+margin)
    new.range[1] <- max(1, new.range[1])
  }

  #Object to return
  R <- list("pheno.file"=pheno.file, "trait.file"= trait.file,
            "range"=range, "trait"=trait, "covariates"=covariates,
            "pheno.transformation"=pheno.transformation,
            "bandwidth"=bandwidth, "smoother"=smoother,
            "z0"=z0, "zmin"=zmin, "chrom"=chrom)


  ####################
  #  Read trait data #
  ####################
  dm <- detect_dm_csv(filename=trait.file, header=TRUE, sep=" ")
  df.laf <- laf_open(dm)
  X <- df.laf[,]
  close(df.laf)
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
  #Permuted trait values
  if(n.perm > 0){
    set.seed(seed)
    perms <- replicate(n=n.perm, expr = {
      sample(x, size=n, replace=FALSE)
    })
  }

  ###################
  # Read pheno data #
  ###################
  #Using the LaF package to read data
  dm <- detect_dm_csv(filename=pheno.file, header=TRUE, sep=" ")
  df.laf <- laf_open(dm)

  #Make sure the trait data is sorted correctly
  if(!all(dm$columns$name[-1] %in% X[[1]])) stop("Not all of the samples in ", pheno.file, " are in ", trait.file, ".\n")
  x <- x[match(dm$columns$name[-1], X[[1]])]
  X <- X[match(dm$columns$name[-1], X[[1]]), ]

  #For saving as we go so we don't use an absurd amount of memory
  if(is.null(temp.prefix)) temp.prefix <- paste0(paste0(sample(c(letters, LETTERS), size=4), collapse=""), "_", chrom)
  tp <- paste0(temp.dir, temp.prefix)
  cat("temp prefix: ", tp, "\n")

  done <- FALSE
  chunk <- 1
  while(!done){
    if(chunk > max(which.chunks)){
      done <- TRUE
      break
    }
    if(!chunk %in% which.chunks){
      chunk <- chunk + 1
      next
    }
    cat("Chunk ", chunk, "\n")
    #read data
    if(chunksize < 0){
      dat <- next_block(df.laf, nrows=nl)
      done <- TRUE
      nstart <- 1
      nend <- nrow(dat)
    }else{
      goto(df.laf, max(1, (chunk-1)*chunksize-margin + 1))
      dat <- next_block(df.laf, nrows=chunksize + 2*margin)
      if(nrow(dat) < chunksize + 2*margin) done <- TRUE
      if(nrow(dat) <= margin) break
      if(chunk==1){
        nstart <- mstart <- 1
      }else{
        nstart <- margin + 1
        mstart <- ceiling(bandwidth/2)
      }
      nmstart <- nstart -mstart + 1
      nend <- min(nrow(dat), nstart + chunksize-1)
      mend <- nrow(dat)-ceiling(bandwidth/2)
      nmend <- nend -mstart  + 1
    }

    if(!is.null(range)){
      if(dat[[1]][nend] < new.range[1]){
        chunk <- chunk + 1
        next
      }
      if(dat[[1]][nend] > new.range[2]){
        done <- TRUE
        nend <- which.min(dat[[1]] < new.range[2])
      }
      if(dat[[1]][nstart] < new.range[1]){
        nstart <- which.max(dat[[1]] >= new.range[1])
      }
    }
    R.temp <- R
    R.temp$chunk <- chunk
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
    R.temp$sts <- stat.func(Y=dat[, -1], x=x,s0=s0)
    R.temp$sts <- data.frame(t(R.temp$sts))
    names(R.temp$sts) <- c("Beta", "SD", "stat")
    R.temp$sts$pos <- dat[[1]]

    if(smoother=="none"){
      R.temp$sts <- R.temp$sts[nstart:nend, ]
      save(R.temp, file=paste0(tp, ".", chunk, ".RData"))
      chunk <- chunk + 1
      next
    }

    ######################################
    # Smooth (non-permutation) statistics#
    ######################################
    cat("Smoothing..\n")
    ys <- smooth.func(x=dat[[1]], y=R.temp$sts$stat,xout=dat[[1]][mstart:mend], bandwidth = bandwidth)
    R.temp$sts.smooth <- data.frame("pos"=dat[[1]][mstart:mend], "ys"=ys)
    R.temp$sts <- R.temp$sts[nstart:nend, ]
    if(is.null(zmin)){
      R.temp$sts.smooth <- R.temp$sts.smooth[nmstart:nmend,]
      save(R.temp, file=paste0(tp, ".", chunk, ".RData"))
      chunk <- chunk + 1
      next
    }
    R.temp$m1 <- mxlist(ys, z0, zmin, pos=dat[[1]][mstart:mend])
    if(!is.null(R.temp$m1)){
      if(s==1) R.temp$m1 <- R.temp$m1[R.temp$m1$mx >= zmin,]
        else R.temp$m1 <- R.temp$m1[R.temp$m1$mx >= zmin[2] & R.temp$m1$mx <= zmin[1],]
      R.temp$m1 <- R.temp$m1[R.temp$m1$pos <= dat[[1]][nend] & R.temp$m1$pos >= dat[[1]][nstart],]
      if(any(R.temp$m1$ix1==1 | R.temp$m1$ix2 ==(mend-mstart + 1))){
        cat("Warning: You may need to increase margin to correctly capture peaks.\n")
      }
      if(nrow(R.temp$m1)==0) R.temp$m1 <- NULL
    }
    R.temp$sts.smooth <- R.temp$sts.smooth[nmstart:nmend,]
    if(n.perm==0){
      save(R.temp, file=paste0(tp, ".", chunk, ".RData"))
      chunk <- chunk + 1
      next
    }

    ###############################
    # Permutation test statistics #
    ###############################
    cat("Calculating permutation statistics...\n")
    N <- nend-nstart + 1
    sum_stat_sq <- rep(0, N)
    sum_stat <- rep(0, N)
    R.temp$mperm <- data.frame(matrix(nrow=0, ncol=6))
    names(R.temp$mperm) <- c("mx", "pos", "start", "stop", "ix1", "ix2")
    for(i in 1:n.perm){
      cat(i, " ")
      sts <- stat.func(dat[,-1], x=perms[,i], s0=s0)
      sm.sts <- smooth.func(x=dat[[1]], y=sts[3,], xout=dat[[1]][mstart:mend], bandwidth = bandwidth)
      sum_stat_sq <- sum_stat_sq + sm.sts[nmstart:nmend]^2
      sum_stat <- sum_stat + sm.sts[nmstart:nmend]
      R.temp$mperm <- rbind(R.temp$mperm, mxlist(sm.sts, z0, zmin, pos=dat[[1]][mstart:mend]))
    }
    cat("\n")
    R.temp$perm.var <- data.frame("pos"=dat[[1]][nstart:nend],
                                  "var"=(sum_stat_sq/(n.perm-1)) - (sum_stat^2)/(n.perm*(n.perm-1)))

    R.temp$mperm <- R.temp$mperm[R.temp$mperm$pos <= dat[[1]][nend] & R.temp$mperm$pos >= dat[[1]][nstart],]
    if(any(R.temp$mperm$ix1==1 | R.temp$mperm$ix2 ==(mend-mstart + 1))){
      cat("Warning: You may need to increase margin to correctly capture peaks.\n")
    }
    cat("Saving..")
    save(R.temp, file=paste0(tp, ".", chunk, ".RData"))
    cat("Saved.\n")
    chunk <- chunk + 1
  }
  rm(dat)
  close(df.laf)

  if(!length(which.chunks)==nchunk){
    return(R)
  }
  cat("Collecting results.\n")
  R <- collect_fret_stats(temp.dir, temp.prefix)
  return(R)
}
