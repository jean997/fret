

#' @title Calculate FRET test statistics
#' @param pheno_file_list Name or vector of names of genomic phenotype files
#' @param trait_file Name of trait file
#' @param mode One of "dry_run", "s0_only", or "full". See description for details.
#' @param s0 Variance inflation constant. If missing, s0 will be estimated
#' using a random chunk of the data. (See description and s0_est_size argument)
#' @param zmin Minimum threshold. If missing, will be set to the 90th percentile of a
#' sample of test statistics (See details below).
#' @param z0 Merging threshold. If missing z0 = 0.3*zmin.
#' @param s0_est_size Number of test statistics used to estimate s0. Alternatively, s0_est_size may be a character
#' string giving a phenotype file name. In this case, all test statistics from that file will be used.
#' @param seed Seed (used for permutations). If missing, a random seed will be chosen and stored.
#' @param n_perm Number of permutations.
#' If n.perm=0, only test statistics for unpermuted data will be calculated.
#' @param pheno_transformation If the phenotype is to be transformed, provide a function taking one argument (the vector of phenotype)
#' and outputting the transformed phenotype.
#' @param trait Name of trait (should match header of trait_file)
#' @param covariates List of covariates to adjust for (names should match header of trait_file)
#' @param stat_type Type of test statistis. May be one of "huber" or "lm".
#' @param huber_maxit Maximum iterations for Huber estimator.
#' @param bandwidth Smoothing bandwidth
#' @param smooother Choice of smoother for smoothing test statistics. Can be
#' one of "ksmooth_0"  or "ksmooth". See details below.
#' @param chunksize Size of chunks to read at one time
#' @param which_chunks Either a vector of chunk numbers or "all". See about running FRET in parallel or on a cluster.
#' @param temp_dir Directory to write temporary chunk output to
#' @param temp_prefix Prefix to use for chunk output
#' @param labels Vector of labels for each phenotype file. These will be used in results tables and also in temporary file names.
#' @param cores Number of cores to use. Using more than one requires the parallel package.
#' @details
#' mode: The function can run in one of three modes. If mode="dry_run" it will report information about the data provided,
#' the model and the number of chunks and then exit. If mode = "s0_only" it will report this information, estimate s0 and zmin and exit.
#' If mode="full" it will proceed to calculate all test statistics and permutation test statistics for the specified chunks.
#'
#' Running on a cluster or in parallel: The which_chunks argument is intended to facilitate breaking a very large job into many small
#' jobs that can be easily submitted to a cluster. It can also help with resuming an analysis that was interrupted. To limit memory
#' requirements, only chunks of size chunksize will be read in and analyzed at one time. Results of these analyses are
#' then written to disk in files named temp_dir/temp_prefix-label.chunknum.RDS. It is important to make sure that temp_dir
#' has enough space to store lots of test statistics. If which_chunks="all", these temporary files will be automatically aggregated
#' into a single set of results. If the analysis is conducted over many jobs, the uster will need to call the collect_fret_stats
#' function to do this themselves(see documentation for collect_fret_stats). In addition to breaking chunks over many nodes
#' or many jobs, the cores parameter can be used to perform calculations using multiple cores via the parallel package.
#'
#' Smoother choice: There are three options for smoothing test statistics. "ksmooth_0" is a box kernel smoother for observations
#' made at integer positions. It assumes that observations at missing positions are equal to 0. This is an appropriate smoother choice
#' for DNase-seq and similar data types. In DNase-seq data, if a position is not present in the data, all samples have 0 cleavages
#' observed at the position so the test statistic is equal to 0. "ksmooth" is a box kernel smoother that assumes observations
#' at missing positions are missing. This is appropriate for bisulfite sequencing data.
#'
#' Estimating s0, zmin, and z0: If s0 is not provided, it will be estimated from the data (if mode = "s0_only" or FALSE).
#' The function will use an amount of data specified by the s0_est_size argument. If this argument is a file name, all the data
#' in that file will be used. If it is an integer, the number of data points specified will be used. If fewer than 1,000,000
#' data points are used, the estimate might be unstable and a warning will be given. If zmin is missing, it will be set to the
#' 90th percentile of the test statistics in the data sample (after correcting using s0). If z0 is missing it will be set
#' to 0.3*zmin.
#'@export
fret_stats <- function(pheno_file_list, trait_file, mode = c("dry_run", "s0_only", "full"),
                       s0, zmin, z0 = if(!missing(zmin)) 0.3*zmin,
                       s0_est_size = pheno_file_list[1], seed, n_perm=0,
                       pheno_transformation=NULL, trait=c("x"), covariates=c(),
                       stat_type=c("huber", "lm"), huber_maxit=50,
                       bandwidth=151, smoother=c("ksmooth_0", "ksmooth"),
                       chunksize=1e5, which_chunks="all",
                       temp_dir ="./",  temp_prefix=NULL, labels=pheno_file_list,
                       cores=1){

  ###########
  # Options #
  ###########
  stopifnot(length(trait)==1)
  if(!missing(pheno_transformation)) stopifnot("function" %in% class(pheno_transformation))

  mode = match.arg(mode)
  stat_type = match.arg(stat_type)
  smoother = match.arg(smoother)
  stopifnot(length(labels)==length(pheno_file_list))
  if(class(s0_est_size)=="character"){
    if(!s0_est_size %in% pheno_file_list) stop(paste0(s0_est_size, " is not in list of phenotype files\n"))
  }else{
    if(!floor(s0_est_size) == s0_est_size) stop(paste0("s0_est_size should be either an integer or a file name. ",
                                                       s0_est_size, "was provided\n"))
  }
  stopifnot(floor(chunksize) == chunksize & chunksize > 0)
  if(class(which_chunks) == "numeric" | class(which_chunks)=="integer"){
    stopifnot(all(floor(which_chunks)==which_chunks))
    stopifnot(all(diff(which_chunks) == 1))
  }else{
    stopifnot(which_chunks=="all")
  }
  #stat type
  use_par <- cores > 1
  if(stat_type=="huber"){
    lm_func <- function(formula, data){
      MASS::rlm(formula, data=data, maxit=huber_maxit)
    }
    if(!use_par) stat_func <- function(Y, x, s0){ huber_stats(Y, x, s0, maxit=huber_maxit)}
      else stat_func <- function(Y, x, s0){ huber_stats_parallel(Y, x, s0=s0, maxit=huber_maxit,
                                                                 cores=cores)}
  }else if(stat_type=="lm"){
    lm_func <- function(formula, data){
      lm(formula, data=data)
    }
    if(!use_par) stat_func <- function(Y, x, s0){ lm_stats_parallel(Y, x, s0=s0, cl=cl)}
    else stat_func <- function(Y, x, s0){ lm_stats_parallel(Y, x, s0=s0, maxit = huber_maxit)}
  }
  #smoother
  stopifnot(floor(bandwidth)==bandwidth & bandwidth > 0)
  smoother <- match.arg(smoother)
  if(smoother=="ksmooth_0"){
    smooth_func <- function(x, y, xout, bandwidth){
      ksmooth_0(x, y,xout, bandwidth)
    }
  }else if(smoother=="ksmooth"){
    smooth_func <- function(x, y, xout, bandwidth){
      ksmooth(x=x, y=y, x.points=xout, bandwidth=bandwidth)$y
    }
  }
  margin <- 5*bandwidth

  if(is.null(temp_prefix)) temp_prefix <- paste0(sample(c(letters, LETTERS), size=4), collapse="")
  if(str_sub(temp_dir, -1) != "/") temp_dir <- paste0(temp_dir, "/")
  cat("temp files will be saved to: ", paste0(temp_dir, temp_prefix, "-label.chunknum.RDS\n"))

  ##############################
  #  Check model, count chunks #
  ##############################

  ######################
  ##  Read trait data ##
  ######################
  dm <- detect_dm_csv(filename=trait_file, header=TRUE, sep=" ")
  df_laf <- laf_open(dm)
  X <- df_laf[,]
  close(df_laf)
  if(!trait %in% names(X)) stop("ERROR: I didn't find colunns matching the specified trait name in the phenotype file.\n")
  if(!all(covariates %in% names(X))) stop("ERROR: I didn't find colunns matching the specified covariates in the phenotype file.\n")
  cat("Trait name: ", trait, "\n")
  #Adjust trait for covariates
  if(length(covariates) > 0){
    cat("Adjusting ", trait, " for ", covariates, "\n")
    ff <- as.formula(paste0(trait, "~", paste0(covariates, collapse="+")))
    fitx <- lm_func(ff, X)
    X[[trait]] <- fitx$residuals
  }
  n <- nrow(X)
  if(n_perm > 0){
    if(missing(seed)) stop("If n_perm > 0, a seed must be provided.\n")
    set.seed(seed)
    perms <- replicate(n=n_perm, expr = {
      sample(1:n, size=n, replace=FALSE)
    })
  }else{
    perms <- NULL
  }
  #########################################
  ## Determine number of chunks per file ##
  #########################################
  chunk_df <- data.frame("File" = pheno_file_list, "label" = labels, "nchunks"=NA,
                         "first_chunk" = NA, "last_chunk" = NA, stringsAsFactors = FALSE)
  chunk_df$first_chunk[1] <- ct <- 1
  for(i in 1:length(pheno_file_list)){
    nl <- determine_nlines(pheno_file_list[i])-1
    chunk_df$nchunks[i] <- max(1, ceiling(nl/chunksize))
  }
  chunk_df$last_chunk <- with(chunk_df, cumsum(nchunks))
  chunk_df$first_chunk <- with(chunk_df, last_chunk-nchunks + 1)
  cat("Chunks per file: \n")
  print(chunk_df)
  total_chunks <- with(chunk_df, sum(nchunks))
  if(class(which_chunks)=="numeric" | class(which_chunks)=="integer"){
    if(any(which_chunks > total_chunks | which_chunks < 1)) stop("ERROR: Some requested chunks are not in 1:", total_chunks, "\n")
    which_files <- with(chunk_df, which( first_chunk <=  max(which_chunks) & last_chunk >= min(which_chunks)))
  }else{
    stopifnot(which_chunks=="all")
    which_chunks <- 1:total_chunks
    which_files <- 1:length(pheno_file_list)
  }
  cat("This job will analyze ", length(which_chunks), " chunks in ", length(which_files), "files.\n")

  if(missing(s0) | missing(zmin)){
    if(!missing(s0)) cat("s0 is given. zmin will be estimated as the 90th percentile of smoothed test statistics in ")
    if(missing(s0) & missing(zmin)) cat("s0 and zmin will be estimated using the data in ")
    if(missing(s0) & !missing(zmin)) warning(paste0("Warning: zmin is provided and s0 is not and will be estimated.",
                                         " Are you sure that's how you want to run?"))
    if(class(s0_est_size)=="character"){
      nc <- with(chunk_df, nchunks[File==s0_est_size])
      cat(s0_est_size, " which contains ", nc, " chunks.\n")
    }else{
      nc <- ceiling(s0_est_size/chunksize)
      if(nc > total_chunks) stop("s0_est_chunk is larger than size of data.\n")
      cat(nc, " random contiguous chunks.\n")
    }
  }


  ### Return if dry run
  if(mode=="dry_run"){
    return(chunk_df)
  }

  #Object to return
  R <- list("pheno_files"=pheno_file_list, "trait_file"= trait_file,
            "trait"=trait, "covariates"=covariates,
            "pheno_transformation"=pheno_transformation, "n_perm"=n_perm,
            "bandwidth"=bandwidth, "smoother"=smoother)

  #####################################
  # Calculate s0, z0, zmin if missing #
  #####################################
  if(missing(s0) | missing(zmin)){
    if(missing(s0)){
      s0 <- 0
      s0miss <- TRUE
    }else{
      s0miss <- FALSE
    }
    cat("Estimating s0 and/or zmin.\n")
    if(class(s0_est_size)=="character"){
      nc <- with(chunk_df, nchunks[File==s0_est_size])
      stats_s0 <- get_stats(s0_est_size, nc, chunksize, margin,
                         stat_func, lm_func, X, trait, covariates, s0=s0,
                         pheno_transformation, "all")
      s0_chunks <- chunk_df %>% filter(File==s0_est_size) %>% with(., first_chunk:last_chunk)
    }else{
      nc <- ceiling(s0_est_size/chunksize)
      start_chunk_s0 <- sample(1:(total_chunks-nc), size=1)
      end_chunk_s0 <- start_chunk_s0 + nc -1
      s0_chunks <- start_chunk_s0:end_chunk_s0
      s0_file_ix <- with(chunk_df, which( first_chunk <=  end_chunk_s0 & last_chunk >= start_chunk_s0))
      stats_s0 <- lapply(s0_file_ix, function(i){
        with(chunk_df, get_stats(File[i], nchunks[i], chunksize, margin,
                  stat_func, lm_func, X, trait, covariates, s0=s0,
                  pheno_transformation,
                  chunks=which(first_chunk[i]:last_chunk[i] %in% s0_chunks)))
      })
      stats_s0 <- do.call(rbind, stats_s0)
    }
    if(s0miss){
      s0 <- with(stats_s0, choose_s0(beta=beta, se=se))
      stats_s0$stat <- with(stats_s0, beta/(se + s0))
      stats_s0$stat[is.na(stats_s0$beta)] <- 0
    }
    if(missing(zmin)){
      stats_s0$stat_smoothed <- with(stats_s0, smooth_func(x=pos, y=stat,xout=pos, bandwidth = bandwidth))
      zmin <- with(stats_s0, as.numeric(quantile(abs(stat_smoothed), 0.9)))
    }
    cat("s0: ", s0, " zmin: ", zmin, "\n")
  }else{
    s0_chunks <- c()
  }
  if(missing(z0)) z0 <- 0.3*zmin

  R$s0 <- s0
  R$zmin <- zmin
  R$z0 <- z0

  if(mode=="s0_only"){
    R$stats_s0 <- stats_s0
    R$s0_chunks <- s0_chunks
    return(R)
  }

  R_temp <- R
  for(chunk in which_chunks){
    file_ix <- with(chunk_df, which( first_chunk <=  chunk & last_chunk >= chunk))
    chunk_in_file <- chunk-chunk_df$first_chunk[file_ix] + 1
    tfile <- paste0(temp_dir, temp_prefix, "-", labels[file_ix], ".", chunk_in_file, ".RDS")
    R_temp$label <- labels[file_ix]
    cat("Analyzing chunk ", chunk, ". Results to be saved to ", tfile, ".\n")
    cat("Calculating test statistics..\n")
    ###############################################
    # Calculate (non-permutation) test statistics #
    ###############################################
    R_temp$stats <- with(chunk_df, get_stats(File[file_ix], nchunks[file_ix], chunksize, margin,
                             stat_func, lm_func, X, trait, covariates, s0=s0,
                             pheno_transformation,
                             chunks=chunk_in_file))
    ######################################
    # Smooth (non-permutation) statistics#
    ######################################
    smooth_ix_range <- c(1, nrow(R_temp$stats))
    peak_pos_range <- R_temp$stats$pos[smooth_ix_range]
    if(!chunk_in_file==1){
      smooth_ix_range[1] <- ceiling(bandwidth/2) + 1
      peak_pos_range[1] <- R_temp$stats$pos[margin + 1]
    }
    if(!chunk_in_file == chunk_df$nchunks[file_ix]){
      smooth_ix_range[2] <- smooth_ix_range[2]-ceiling(bandwidth/2)
      peak_pos_range[2] <- with(R_temp$stats, pos[ length(pos)-margin -1])
    }
    cat("Smoothing..\n")
    R_temp$stats$stat_smoothed <- with(R_temp$stats, smooth_func(x=pos, y=stat,xout=pos, bandwidth = bandwidth))
    R_temp$stats <- R_temp$stats[ smooth_ix_range[1]:smooth_ix_range[2],]
    R_temp$peaks <- with(R_temp$stats, mxlist(stat_smoothed, z0, zmin, pos=pos)) %>%
                   filter(mx > zmin) %>%
                   filter(pos >= peak_pos_range[1] & pos <= peak_pos_range[2])

    if(n_perm==0){
      saveRDS(R_temp, file=tfile)
      next
    }
    ###############################
    # Permutation test statistics #
    ###############################
    perm <-with(chunk_df, fret:::get_perm_peaks(File[file_ix], chunk_in_file, perms,
                                        nchunks[file_ix], chunksize, margin, smooth_ix_range, peak_pos_range,
                                        stat_func, lm_func, smooth_func,  bandwidth,
                                        X, trait, covariates, s0, z0, zmin,
                                        pheno_transformation=NULL))
    R_temp$perm_peaks <- perm$perm_peaks
    R_temp$stats$perm_var <- perm$perm_var[smooth_ix_range[1]:smooth_ix_range[2]]
    R_temp$range <- peak_pos_range
    saveRDS(R_temp, file=tfile)
  }
  #if(!length(which.chunks)==nchunk){
    return(R)
  #}
  #cat("Collecting results.\n")
  #R <- collect_fret_stats(temp_dir, temp_prefix, which.chunk = which_chunks)
  #return(R)
}
