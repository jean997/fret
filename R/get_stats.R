#'@title Helper function for calculating test statistics
#'@param file_name Name of phenotyp file
#'@param nchunks Number of chunks in file
#'@param chunksize Chunk size
#'@param margin Margin
#'@param X Data frame of trait information
#'@param trait String; name of trait in X
#'@param covariates Vector of names of covariates in X
#'@param sample String; name of sample name in X
#'@param s0 s0
#'@param stat_fun Function to calculate test stats
#'@param resid_fun Function to calculate residual
#'@param cores
#'@param pheno_transformation Phenotype transformation
#'@param chunks Which chunks; may only be "all" or contiguous integers
#'@export
get_stats <- function(file_name, nchunks, chunksize, margin,
                       X, trait, covariates, sample, s0, stat_fun, resid_fun,
                      cores, libs,
                      pheno_transformation=NULL, chunks="all"){
  #Check chunk argument
  if(class(chunks) == "numeric" | class(chunks)=="integer"){
    stopifnot(all(chunks %in% 1:nchunks))
  }else{
    stopifnot(chunks=="all")
    chunks <- 1:nchunks
  }
  chunks <- sort(chunks)
  stopifnot(all(diff(chunks)==1)) #only contiguous chunks

  #Open pheno file
  dm <- detect_dm_csv(filename=file_name, header=TRUE, sep=" ")
  df_laf <- laf_open(dm)

  #Make sure the trait data is sorted correctly
  if(!all(dm$columns$name[-1] %in% X[[sample]])){
    nms <- dm$columns$name[-1]
    cat(paste0(nms[!nms %in% X[[sample]]], collapse=" "), "\n")
    cat(paste0(X[[sample]], collapse=" "), "\n")
    stop("Not all of the samples in ", file_name, " are in the trait file.\n")
  }
  X <- X[match(dm$columns$name[-1], X[[sample]]),]

  lmargin <- rep(0, length(chunks))
  rmargin <- rep(0, length(chunks))
  if(!min(chunks)==1) lmargin[1] <- margin
  if(!max(chunks)==nchunks) rmargin[length(chunks)] <- margin

  goto(df_laf, max(1, (chunks[1]-1)*chunksize-lmargin[1] + 1))
  stats <- lapply(1:length(chunks), function(i){
      dat <- next_block(df_laf, nrows=chunksize + lmargin[i] + rmargin[i])
      sts <- calc_stats(dat, X, s0, trait, covariates, pheno_transformation, stat_fun, resid_fun, cores, libs)
      strt <- max(1, (chunks[i]-1)*chunksize-lmargin[i] + 1)
      ix <- strt:(strt + nrow(sts)-1)
      cbind(sts, ix, dat[[1]])
  })
  stats <- as.data.frame(do.call(rbind, stats))
  names(stats) <- c("beta", "se", "stat", "index", "pos")
  close(df_laf)
  return(stats)
}

calc_stats <- function(dat, X, s0,
                       trait, covariates, pheno_transformation,
                       stat_fun, resid_fun, cores, libs){

  if(!is.null(pheno_transformation)){
    cat("Adjusting phenotype.\n")
    dat[,-1] <- apply(dat[,-1], MARGIN=2, FUN=function(y){
      pheno_transformation(y)
    })
    if(any(is.na(dat[,-1])) | any(!is.finite(unlist(dat[,-1])))){
      stop("ERROR: Phenotype adjustment gives infinite or missing value.\n")
    }
  }
  if(length(covariates) > 0){
    cat("Regressing covariates on phenotype.\n")
    dat[,-1] <- apply(dat[,-1], MARGIN=1, FUN=function(y){
      resid_fun(y, X, covariates)
    })
  }
  cat("Calculating stats..\n")
  stats <- stats1(Y=dat[,-1], x=X[[trait]], s0,  cores, stat_fun, libs)
  return(stats)
}

