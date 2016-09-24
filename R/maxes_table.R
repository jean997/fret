#'@export
fret_stats <- function(dat.file, pheno.file, s0, seed, n.perm, zmin=NULL,
                        pheno.transformation=NULL, trait=c("x"), covariates=c(),
                        z0=zmin*0.3, save.perm.stats=FALSE, range=NULL, chunksize=10000,
                        bandwidth=50, buffer=2*bandwidth, out.file=NULL, chrom="chr1",
                        stat.type=c("huber", "lm"), maxit=50, smoother=c("ksmooth_0", "ksmooth")){
  #Options
  if(!length(trait) ==1) stop("ERROR: Handling of multivariate traits is not implemented yet!\n")
  if(is.null(zmin)) stopifnot(n.perm==0)
  #stat type
  stat.type <- match.arg(stat.type)
  if(stat.type=="huber"){
    lm.func <- function(formula, data){
      rlm(formula, data=data, maxit=maxit)
    }
    stat.func <- function(Y, x, s0){ huber_stats(Y, x, s0, maxit=maxit)}
  }else if(staty.type=="lm"){
    lm.func <- function(formula, data){
      lm(formula, data=data)
    }
    stat.func <- lm_stats
  }
  #smoother type
  smoother <- match.arg(smoother)
  if(smoother=="ksmoooth_0"){
    smooth.func <- function(x, y, bandwidth){ksmooth_0_old(x, y,x, bandwidth)}
  }else{
    smooth.func <- function(x, y, bandwidth){ ksmooth(x=x, y=y, x.points=x, bandwidth=bandwidth)$y}
  }

  stopifnot(length(z0)==length(zmin) | is.null(zmin))
  s <- length(zmin)
  stopifnot(s %in% c(0, 1, 2))
  if(!is.null(pheno.transformation)) stopifnot("function" %in% class(pheno.transformation))
  chunksize <- floor(chunksize)
  #Read phenotype
  X <- read_delim(pheno.file, delim=" ")
  if(!all(trait %in% names(X))) stop("ERROR: I didn't find colunns matching the specified trait name in the phenotype file.\n")
  if(!all(covariates %in% names(X))) stop("ERROR: I didn't find colunns matching the specified covariates in the phenotype file.\n")
  n <- nrow(X)
  #Adjust phenotype for covariates
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

  #Read data
  #Read header
  h <- readLines(dat.file, n=1)
  h <- unlist(strsplit(h, " "))
  if(! h[1]=="pos") stop("First column of data matrix should be genomic position and be named 'pos'")

  #Make sure the phenotype is sorted correctly
  if(!all(h[-1] %in% X[[1]])) stop("Not all of the samples in the data file are in the phenotype file.\n")
  x <- x[match(h[-1], X[[1]])]
  X <- X[match(h[-1], X[[1]]), ]
  #Range
  if(!is.null(range)){
    new.range <- c(range[1]-ceiling(bandwidth/2)-buffer , range[2]+ceiling(bandwidth/2)+buffer)
    new.range[1] <- max(1, new.range[1])
    cat("Expanded range: ", new.range, "\n")
    dat <- read_data_range1(dat.file, new.range, chunksize=chunksize)
    pos <- dat$pos
    ix1 <- min(which(pos >=range[1]))
    ix2 <- max(which(pos <= range[2]))
  }else{
    col_types=paste0(c("i", rep("d", length(h)-1)), collapse="")
    dat <- read_delim(dat.file, col_types=col_types, delim=" ")
    pos <- dat$pos
    ix1 <- 1
    ix2 <- nrow(dat)
  }
  pos.out = pos[ix1:ix2]
  cat(dim(dat)[1], " positions.\n")
  #Adjust phenotype for covariates
  if(!is.null(pheno.transformation)){
    cat("Adjusting phenotype.\n")
    yadj <- apply(dat[,-1], MARGIN=2, FUN=function(y){
      pheno.transformation(y)
    })
    if(any(is.na(yadj)) | any(!is.finite(yadj))) stop("ERROR: Phenotype adjustment gives infinite or missing value.\n")
    dat[,-1] <- yadj
  }
  if(length(covariates) > 0){
    cat("Regressing covariates on phenotype.\n")
    yadj <- apply(dat[,-1], MARGIN=1, FUN=function(y){
      ff <- as.formula(paste0("y~", paste0(covariates, collapse="+")))
      fity <- lm.func(ff, X)
      fity$residuals
    })
    dat[,-1] <- yadj
  }
  #Calculate statistics
  cat("Calculating test statistics..\n")
  sts <- stat.func(Y=dat[, -1], x=x,s0=s0)
  sts <- data.frame(t(sts))
  names(sts) <- c("Beta", "SD", "stat")
  sts$pos <- pos
  #Smooth statistics
  ys <- smooth.func(x=pos, y=sts$stat, bandwidth = bandwidth)
  stats.smooth <- data.frame("pos"=pos.out, "ys"=ys[ix1:ix2])

  if(is.null(zmin)){
    R <- list("dat.file"=dat.file, "pheno.file"= pheno.file,
              "range"=range, "trait"=trait, "covariates"=covariates,
              "pheno.transformation"=pheno.transformation,
              "stats" = sts, "stats.smooth"=stats.smooth)
    if(!is.null(out.file)){
      save(R, file=out.file)
      return(NULL)
    }
    return(R)
  }

  #Find peak heights
  #Intervals defined by z0:
  max1 <- mxlist(ys, z0, zmin, return.ix = TRUE)
  if(s==1) max1 <- max1[max1$mx > zmin,]
    else max1 <- max1[max1$mx > zmin[1] | max1$mx < zmin[2],]
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
    R <- list("max1"=max1,"file"=dat.file,
              "range"=range, "trait"=trait, "covariates"=covariates,
              "pheno.transformation"=pheno.transformation,
              "stats" = sts, "stats.smooth"=stats.smooth,
              "z0"=z0, "zmin"=zmin)
    if(!is.null(out.file)){
      save(R, file=out.file)
      return(NULL)
    }
    return(R)
  }

  #Permutation stats
  cat("Calculating permutation statistics...\n")
  max.perm <- data.frame("mx"=c(), "ix"=c(), "iv1"=c(), "iv2"=c())
  sts.perm <- apply(perms, MARGIN=2, FUN=function(l){
   stat.func(dat[,-1], x=l, s0=s0)
  })
  stats.perm <- sts.perm[ seq(3, nrow(sts.perm), by=3),]
  beta.perm <- sts.perm[ seq(1, nrow(sts.perm), by=3),]
  sd.perm <- sts.perm[ seq(2, nrow(sts.perm), by=3),]
  stats.perm.smooth <- apply(stats.perm, MARGIN=2, FUN=function(yy){
    smooth.func(x=pos, y=yy, bandwidth = bandwidth)
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
    if(!is.null(range)) max.perm <- max.perm[max.perm$mx >= range[1] & max.perm$mx <= range[2],]
  }
  if(save.perm.stats){
    R <- list("max1"=max1, "max.perm"=max.perm, "file"=dat.file,
              "z0"=z0, "zmin"=zmin, "n.perm"=n.perm, "perm.var"=vv,
              "range"=range, "trait"=trait, "covariates"=covariates,
              "pheno.transformation"=pheno.transformation,
              "stats"=sts, "stats.smooth"=stats.smooth,
              "stats.perm" = cbind(pos, stats.perm),
              "beta.perm" = cbind(pos, beta.perm),
              "sd.perm" = cbind(pos, sd.perm),
              "stats.perm.smooth"= cbind(pos, stats.perm.smooth))
  }else{
    R <- list("max1"=max1, "max.perm"=max.perm, "file"=dat.file,
              "z0"=z0, "zmin"=zmin, "n.perm"=n.perm, "perm.var"=vv,
              "range"=range, "trait"=trait, "covariates"=covariates,
              "pheno.transformation"=pheno.transformation,
              "stats"=sts, "stats.smooth" = stats.smooth)
  }
  if(!is.null(out.file)){
    save(R, file=out.file)
    return(nrow(max.perm))
  }
  return(R)
}

read_data_range0 <- function(dat.file, range, chunksize=10000){
  con <- file(dat.file, "rb")
  h <- read_lines(con, n_max=1)
  h <- unlist(strsplit(h, " "))
  pos <- 0
  while(pos < range[1]){
    dat <- read_delim(con, n_max=chunksize)
    pos <- dat[chunksize, 1]
  }
  while(pos < range[2]){
    dd <- read_delim(con, n_max=chunksize)
    pos <- dat[chunksize, 1]
    dat <- rbind(dat, dd)
  }
  names(dat) <- h
  ix1 <- min(which(dat$pos >= range[1]))
  ix2 <- max(which(dat$pos <= range[2]))
  return(dat[ix1:ix2,])
}

#Use this until readr gets fixed?
#Right now readr always dumps you at the end even if you only read part of the file
read_data_range1 <- function(dat.file, range, chunksize=10000){
  con <- file(dat.file, "rb")
  h <- readLines(con, n=1)
  h <- unlist(strsplit(h, " "))
  top <- seek(con)
  pos <- 0
  ct <- 1
  col_types=paste0(c("i", rep("d", length(h)-1)), collapse="")
  while(pos < range[1]){
    cat(ct, " ")
    dat <- read_delim(con, n_max=chunksize, skip=(ct-1)*chunksize,
                      col_names = h, col_types=col_types,
                      delim=" ")
    if(nrow(dat)==0) stop("ERROR: Range begins after end of file.\n")
    pos <- dat[[1]][nrow(dat)]
    cat(pos, ".. ")
    seek(con, top)
    ct <- ct + 1
  }
  while(pos < range[2]){
    cat(ct, " ")
    dd <- read_delim(con, n_max=chunksize, skip=(ct-1)*chunksize,
                     col_names = h, col_types=col_types,
                     delim=" ")
    if(nrow(dd)==0){
      pos <- Inf
    }else{
      pos <- dd[[1]][nrow(dd)]
      dat <- rbind(dat, dd)
    }
    cat(pos, ".. ")
    ct <- ct + 1
    seek(con, top)
  }

  ix1 <- min(which(dat$pos >= range[1]))
  ix2 <- max(which(dat$pos <= range[2]))
  close(con)
  return(dat[ix1:ix2,])
}
