#'@export
get_perm_peaks <- function(file_name, chunk, perms,
                           nchunks, chunksize, margin, smooth_ix_range, peak_pos_range,
                           stat_fun, resid_fun, smooth_func, cores, libs, bandwidth,
                           X, trait, covariates, s0, z0, zmin,
                           pheno_transformation=NULL){

  #Check chunk argument
  stopifnot(class(chunk)=="numeric" | class(chunk)=="integer")

  #Open pheno file
  dm <- detect_dm_csv(filename=file_name, header=TRUE, sep=" ")
  df_laf <- laf_open(dm)

  #Make sure the trait data is sorted correctly
  if(!all(dm$columns$name[-1] %in% X[[1]])){
    stop("Not all of the samples in ", file_name, " are in the trait file.\n")
  }
  X <- X[match(dm$columns$name[-1], X[[1]]),]

  if(chunk==1) lmargin <- 0
    else lmargin <- margin
  if(chunk==nchunks) rmargin <- 0
    else rmargin <- margin

  goto(df_laf, max(1, (chunk-1)*chunksize-lmargin + 1))
  dat <- next_block(df_laf, nrows=chunksize + lmargin + rmargin)
  close(df_laf)
  #If necessary, transform phenotype
  if(!is.null(pheno_transformation)){
    cat("Adjusting phenotype.\n")
    dat[,-1] <- apply(dat[,-1], MARGIN=2, FUN=function(y){
      pheno_transformation(y)
    })
    if(any(is.na(dat[,-1])) | any(!is.finite(dat[,-1]))) stop("ERROR: Phenotype adjustment gives infinite or missing value.\n")
  }
  if(length(covariates) > 0){
    cat("Regressing covariates on phenotype.\n")
    dat[,-1] <- apply(dat[,-1], MARGIN=1, FUN=function(y){
      resid_fun(y, X, covariates)
    })
  }
  pos_out <- dat$pos[smooth_ix_range[1]:smooth_ix_range[2]]
  cat("Calculating test statistics for ", ncol(perms), " permutations.\n")

  perm_traits <- apply(perms, 2, function(ix){X[[trait]][ix]})
  perm_stats <- stats_many(Y=dat[,-1], X = perm_traits, s0, cores, stat_fun, libs)
  perm_stats_sm <- apply(perm_stats, 2, function(sts){
    smooth_func(x=dat$pos, y=sts, xout=pos_out, bandwidth=bandwidth)
  })
  if(!sum(is.na(perm_stats_sm))==0){
    cols <- unique(which(is.na(perm_stats_sm), arr.ind = TRUE)[,2])
    cat("NAs found: ", chunk, cols, "\n")
  }
  cat("Calculating variance of smoothed permutation statistics.\n")
  perm_var <- apply(perm_stats_sm[pos_out >= peak_pos_range[1] &  pos_out <= peak_pos_range[2], ], 1, var, na.rm=TRUE)
  cat("Finding peaks in smoothed permutation statistics.\n")
  perm_peaks <- apply(perm_stats_sm, 2, function(sm_sts){
    mxlist(sm_sts, z0, zmin, pos=pos_out)
  })
  perm_peaks <- as.data.frame(do.call(rbind, perm_peaks))
  names(perm_peaks) <- c("mx", "pos", "start", "stop", "ix1", "ix2")
  close(df_laf)
  return(list("perm_var"=perm_var, "perm_peaks"=perm_peaks))
}
