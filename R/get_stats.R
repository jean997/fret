#'@export
get_stats <- function(file_name, nchunks, chunksize, margin,
                      stat_func, lm_func, X, trait, covariates, s0,
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
  if(!all(dm$columns$name[-1] %in% X[[1]])){
    stop("Not all of the samples in ", file_name, " are in the trait file.\n")
  }
  X <- X[match(dm$columns$name[-1], X[[1]]),]

  lmargin <- rep(0, length(chunks))
  rmargin <- rep(0, length(chunks))
  if(!min(chunks)==1) lmargin[1] <- margin
  if(!max(chunks)==nchunks) rmargin[length(chunks)] <- margin

  goto(df_laf, max(1, (chunks[1]-1)*chunksize-lmargin[1] + 1))
  stats <- lapply(1:length(chunks), function(i){
      dat <- next_block(df_laf, nrows=chunksize + lmargin[i] + rmargin[i])
      sts <- calc_stats(dat, X,trait, covariates, pheno_transformation, stat_func, lm_func, s0)
      strt <- max(1, (chunks[i]-1)*chunksize-lmargin[i] + 1)
      ix <- strt:(strt + nrow(sts)-1)
      cbind(sts, ix, dat[[1]])
  })
  stats <- as.data.frame(do.call(rbind, stats))
  names(stats) <- c("beta", "se", "stat", "index", "pos")
  close(df_laf)
  return(stats)
}

calc_stats <- function(dat, X, trait, covariates, pheno_transformation, stat_func, lm_func, s0){

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
      ff <- as.formula(paste0("y~", paste0(covariates, collapse="+")))
      suppressWarnings(fity <- lm_func(ff, X))
      fity$residuals
    })
  }
  cat("Calculating stats..\n")
  stats <- t(stat_func(Y=dat[, -1], x=X[[trait]],s0=s0))
  return(stats)
}

