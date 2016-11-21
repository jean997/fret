
#' Calculate FRET test statistics
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
fret_stats2 <- function(pheno.file, trait.file, s0, seed, n.perm, zmin=NULL, z0=zmin*0.3,
                        pheno.transformation=NULL, trait=c("x"), covariates=c(),
                        save.perm.stats=FALSE, range=NULL, chunksize=1e5,
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
  con <- file(pheno.file, "rb")
  #Read header
  h <- readLines(con, n=1)
  h <- unlist(strsplit(h, " "))
  #if(! h[1]=="pos") stop("First column of data matrix should be genomic position and be named 'pos'")

  #Make sure the trait data is sorted correctly
  if(!all(h[-1] %in% X[[1]])) stop("Not all of the samples in ", pheno.file, " are in ", trait.file, ".\n")
  x <- x[match(h[-1], X[[1]])]
  X <- X[match(h[-1], X[[1]]), ]

  #For now read in the entire file
  col_types=paste0(c("i", rep("d", length(h)-1)), collapse="")
  dat <- read_delim(con, col_types=col_types, delim=" ", col_names=header)
  close(con)
  nn <- nrow(dat)
  #Range
  if(!is.null(range)){
    #We need to look at data beyond the specified range to get the smoothing right
    #And also because some peaks might go outside of the range
    new.range <- c(range[1]-3*bandwidth , range[2]+3*bandwidth)
    new.range[1] <- max(1, new.range[1])
    if(dat[[1]][nn] < new.range[1]) stop("ERROR: Range begins after end of file.\n")
    ix1 <- which.max(dat[[1]] >= new.range[1])
    ix2 <- which.min(dat[[1]] <= new.range[2])
    dat <- dat[ix1:ix2,]
    nn <- nrow(dat)
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
  sts <- stat.func(Y=dat[, -1], x=x,s0=s0)
  sts <- data.frame(t(my.sts))
  names(sts) <- c("Beta", "SD", "stat")
  sts$pos <- dat[[1]]

  if(smoother=="none" & n.perm==0){
    R <- list("pheno.file"=pheno.file, "trait.file"= trait.file,
              "range"=range, "trait"=trait, "covariates"=covariates,
              "pheno.transformation"=pheno.transformation,
              "bandwidth"=bandwidth, "smoother"=smoother,
              "stats" = sts)
    if(!is.null(out.file)){
      save(R, file=out.file)
      return(NULL)
    }
    return(R)
  }

  ######################################
  # Smooth (non-permutation) statistics#
  ######################################
  if(!is.null(range)){
    ix1 <- which.max(sts$pos >= range[1])
    ix2 <- which.min(sts$pos <= range[2])
    pos.out <- sts$pos[ix1:ix2]
  }else{
    pos.out <- sts$pos
  }
  ys <- smooth.func(x=pos, y=sts$stat,xout=pos.out, bandwidth = bandwidth)
  sts.smooth <- data.frame("pos"=pos.out, "ys"=ys[ix1:ix2])

  if(is.null(zmin)){
    R <- list("pheno.file"=pheno.file, "trait.file"= trait.file,
              "range"=range, "trait"=trait, "covariates"=covariates,
              "pheno.transformation"=pheno.transformation,
              "bandwidth"=bandwidth, "smoother"=smoother,
              "stats" = sts, "stats.smooth"=sts.smooth)
    if(!is.null(out.file)){
      save(R, file=out.file)
      return(NULL)
    }
    return(R)
  }

  #######################################
  # Find (non-permutation) peak heights #
  #######################################

  max1 <- mxlist(ys, z0, zmin, return.ix = TRUE)
  if(s==1) max1 <- max1[max1$mx > zmin,] #Unsigned
  else max1 <- max1[max1$mx > zmin[1] | max1$mx < zmin[2],] #signed
  if(nrow(max1) == 0){
    cat("No peaks exceed ", zmin, "\n")
  }
  #Convert index to position
  max1$pos <- pos[max1$ix]
  max1$iv1 <- pos[max1$iv1]
  max1$iv2 <- pos[max1$iv2]
  max1$chr <- chrom
  max1 <- max1[, c("mx", "chr", "pos", "iv1", "iv2")]
  if(!is.null(range)) max1 <- max1[max1$mx >= range[1] & max1$mx <= range[2],]
  if(n.perm==0){
    R <- list("pheno.file"=pheno.file, "trait.file"= trait.file,
              "range"=range, "trait"=trait, "covariates"=covariates,
              "pheno.transformation"=pheno.transformation,
              "bandwidth"=bandwidth, "smoother"=smoother,
              "z0"=z0, "zmin"=zmin,
              "stats" = sts, "stats.smooth"=sts.smooth, "max1"=max1)
    if(!is.null(out.file)){
      save(R, file=out.file)
      return(NULL)
    }
    return(R)
  }



  #########################################
  # Calculate permutation test statistics #
  #########################################
  #For saving as we go so we don't use an absurd amount of memory
  tdir <- paste0(tmp.dir, paste0(sample(c(letters, LETTERS), size=4), collapse=""), "_", chrom)
  system(paste0("mkdir ", tdir))

  cat("Calculating permutation statistics...\n")
  max.perm <- data.frame("mx"=c(), "ix"=c(), "iv1"=c(), "iv2"=c())
  for(i in 1:n.perm){
    cat(i, " ")
    #stats
    sts.perm <- stat.func(Y=dat[, -1], x=perms[,i],s0=s0)
    #smooth
    ys <- smooth.func(x=pos, y=sts.perm$stat,xout=pos.out, bandwidth = bandwidth)
    sts.perm.smooth <- data.frame("pos"=pos.out, "ys"=ys[ix1:ix2])
    #save
    if(save.perm.stats){
      save(sts.perm, file=paste0(tdir, "/stats_", i, ".RData"))
    }
    save(sts.perm.smooth, file=paste0(tdir, "/stats_sm_", i, ".RData"))
    #peaks
    maxperm <- mxlist(ys, z0, zmin, return.ix = TRUE)
    if(s==1) maxperm <- maxperm[maxperm$mx > zmin,] #Unsigned
      else maxperm <- maxperm[maxperm$mx > zmin[1] | maxperm$mx < zmin[2],] #signed
    max.perm <- rbind(max.perm, maxperm)
  }

  if(nrow(max.perm) > 0){
    #Convert positions
    max.perm$pos <- pos.out[max.perm$ix]
    max.perm$iv1 <- pos.out[max.perm$iv1]
    max.perm$iv2 <- pos.out[max.perm$iv2]
    max.perm$chr <- chrom
    max.perm <- max.perm[order(max.perm$pos, decreasing = FALSE), c("mx", "chr", "pos", "iv1", "iv2")]
    if(!is.null(range)) max.perm <- max.perm[max.perm$pos >= range[1] & max.perm$pos <= range[2],]
  }

  R <- list("pheno.file"=pheno.file, "trait.file"= trait.file,
            "range"=range, "trait"=trait, "covariates"=covariates,
            "pheno.transformation"=pheno.transformation,
            "bandwidth"=bandwidth, "smoother"=smoother,
            "z0"=z0, "zmin"=zmin,
            "stats" = sts, "stats.smooth"=sts.smooth,
            "max1"=max1, "max.perm"=max.perm)
  save(R, file=paste0(tdir, "/R.RData"))
  ###########################################
  # Calculate variance of permutation stats #
  ###########################################
  rm(dat, sts, sts.smooth)
  sts.perm.smooth <- getobj(paste0(tdir, "/stats_sm_", 1, ".RData"))
  unlink(paste0(tdir, "/stats_sm_", 1, ".RData"))
  for(i in 2:n.perm){
    ss <- getobj(paste0(tdir, "/stats_sm_", i, ".RData"))
    sts.perm.smooth <- cbind(sts.perm.smooth, ss)
    unlink(paste0(tdir, "/stats_sm_", i, ".RData"))
  }
  vv <- apply(sts.perm.smooth, MARGIN=1, FUN=var)

  ########
  # Save #
  ########
  if(save.perm.stats){
    sts.perm <- getobj(paste0(tdir, "/stats_", 1, ".RData"))
    unlink(paste0(tdir, "/stats_", 1, ".RData"))
    for(i in 2:n.perm){
      ss <- getobj(paste0(tdir, "/stats_", i, ".RData"))
      sts.perm <- cbind(sts.perm, ss)
      unlink(paste0(tdir, "/stats_", i, ".RData"))
    }
    R <- getobj(paste0(tdir, "/R.RData"))
    unlink(paste0(tdir, "/R.RData"))
    R$stats.perm <- sts.perm
    R$stats.perm.smooth <- sts.perm.smooth
    R$perm.var <- vv
  }else{
    R <- getobj(paste0(tdir, "/R.RData"))
    unlink(paste0(tdir, "/R.RData"))
    R$perm.var=vv
  }
  if(!is.null(out.file)){
    save(R, file=out.file)
    return(nrow(max.perm))
  }
  return(R)
}
