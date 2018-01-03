#' Calculate per-base error rates for original and permutation peaks
#'@param fret_obj fret object or a file name or a vector of file names.
#'@param segment_lengths data frame with (at least) two columns: segment and length in base-pairs
#'@param segment_bounds data frame with (at least) two columns: start and stop base pair position
#'@export
fret_rates_prelim <- function(fret_obj, segment_lengths, segment_bounds){


  if(class(fret_obj)=="character"){
    dat <- lapply(fret_obj, function(file){
        f <- readRDS(file)
        f$peaks$label <- f$label
        list("peaks" = f$peaks, "perm_peaks" = f$perm_peaks,
             "range" = f$range, "n_perm" = f$n_perm, "zmin" = f$zmin)
    })
    peaks <- do.call(rbind, lapply(dat, function(x){x$peaks}))
    perm_peaks <- do.call(rbind, lapply(dat, function(x){x$perm_peaks}))
    size <- do.call(sum, lapply(dat, function(x){x$range[2]-x$range[1] + 1}))
    n_perm <- sapply(dat, function(x){x$n_perm})
    stopifnot(all(n_perm==n_perm[1]))
    n_perm <- n_perm[1]
    zmin <- sapply(dat, function(x){x$zmin})
    stopifnot(all(zmin==zmin[1]))
    zmin <- zmin[1]
  }else{
    peaks <- fret_obj$peaks
    peaks$label <- fret_obj$label
    perm_peaks <- fret_obj$perm_peaks
    size <- fret_obj$range[2]-fret_obj$range[1] + 1
    n_perm <- fret_obj$n_perm
    zmin <- fret_obj$zmin
  }

  if(!missing(segment_bounds)){
    if(!missing(segment_lengths)) warning("both segment_bounds and segment_lengths were provided. segment_lengths will be ignored.\n")
    if(!all(segment_bounds$chrom == segment_bounds$chrom[1])) stop("Please only provide one chromosome at a time.\n")
    k <- nrow(segment_bounds)
    gaps <- segment_bounds$start[-1]-segment_bounds$stop[-k]
    stopifnot(all(gaps >=0))
    stopifnot(all(segment_bounds$stop > segment_bounds$start))
    brks <- sort(unlist(segment_bounds[, c("start", "stop")]))

    segs <- findInterval(peaks$pos, brks)
    if(!all(segs %%2 == 1)){
      ix <- which(segs%%2==0)
      if(all(perm_peaks$pos[ix] %in% segment_bounds$stop)){
        segs[ix] <- segs[ix]-1
      }else{
        stop("Some permutation peaks are not in given segments.\n")
      }
    }
    peaks$segment <- (segs+1)/2

    segs <- findInterval(perm_peaks$pos, brks)
    if(!all(segs %%2 == 1)){
      ix <- which(segs%%2==0)
      if(all(perm_peaks$pos[ix] %in% segment_bounds$stop)){
        segs[ix] <- segs[ix]-1
      }else{
        stop("Some permutation peaks are not in given segments.\n")
      }
    }
    perm_peaks$segment <- (segs+1)/2
    segment_bounds$length <- with(segment_bounds, stop-start+1)
    segment_bounds$segment <- seq(nrow(segment_bounds))
    segment_lengths <- segment_bounds
    cat(segment_bounds$chrom[1], "\n")
  }


  has_segment_main <- "segment" %in% names(peaks)
  has_segment_peak <- "segment" %in% names(perm_peaks)
  stopifnot(has_segment_main==has_segment_peak)
  if(!has_segment_main){
    peaks$segment <- 1
    perm_peaks$segment <- 1
    if(missing(segment_lengths)){
      segment_lengths <- data.frame("segment" = 1, length=size)
    }
  }
  stopifnot(all(c("segment", "length") %in% names(segment_lengths)))

  npeaks <- peaks %>% group_by(segment) %>% summarize(np=n())
  nperm_peaks <- perm_peaks %>% group_by(segment) %>% summarize(npp=n())
  n <- full_join(npeaks, nperm_peaks, by="segment") %>%
    full_join(., segment_lengths, by="segment") %>%
  mutate("np" = recode(np, .missing=as.integer(0)),
         "npp" = recode(npp, .missing=as.integer(0)))
  if(any(is.na(n$length))) stop("Not all segment lengths were provided.\n")


  n$max_lambda_perbase <- 0
  #For each peak in peaks and perm_peaks we want to find the value of lambda
  #at the max height of the peak.
  #We also want to find the overall maximum achievable lambda in each segement.
  #This is important since some segments have a very low max lambda
  #(e.g. we almost never see a permutation peak above z0 in these segments)

  peaks$lambda_perbase <- 0

  for(i in seq(n$segment)){
    k <- n$segment[i]
    cat(k, "..")
    #No or almost no permutation peaks above z0 --> max_lambda = 0; lambda per base of any peak in data is 0
    if(n$npp[i] <2 ){
      n$max_lambda_perbase[i] = 0
      peaks <- mutate(peaks, lambda_perbase = if_else(segment==k, true=0, false=lambda_perbase))
      next
    }

    lpb <-  filter(perm_peaks, segment==k) %>% dplyr::select(mx) %>%
            mutate( "mx" = sort(mx, decreasing=TRUE),
              "lambda_per_base" = (seq(n$npp[i])/(n_perm*n$length[i])))

    #Find maximum possible lambda value (i.e. rate at zmin)
    n$max_lambda_perbase[i] <- get_rate_with_thresh(zmin, lpb, num_points_project=10)

    if(n$np[i] > 0){
      rts <- filter(peaks, segment==k) %>% dplyr::select(mx) %>%
             unlist(.) %>%
             get_rate_with_thresh(thresh = ., lpb, num_points_project=10)
      peaks$lambda_perbase[peaks$segment==k] <- rts
    }
  }
  return(list("peaks" = peaks, "segment_info" = n))
}


#If thresh is larger than the largest permutation peak then we extrapolate lambda per base using
# the top num_points_project points and linear regression
# Otherwise we use loess

get_rate_with_thresh <- function(thresh, lpb, num_points_project=10){
  stopifnot(all(diff(lpb$mx) <= 0))
  stopifnot(all(diff(lpb$lambda_per_base) >= 0))
  if(nrow(lpb) < num_points_project){
    warning("Fewer permutations than desired for projection. Consider running more permutations?\n")
    num_points_project <- nrow(lpb)
  }
  #fit_lo <- loess(log10(lambda_per_base) ~ mx, data=lpb)
  fit_sp <- with(lpb, splinefun(x=mx, y=log10(lambda_per_base)))
  preds <- 10^(fit_sp(thresh))
  #if(any(is.na(preds) & thresh < max(lpb$mx))) stop("Something went wrong with loess smoothing..\n")
  if(any(thresh > lpb$mx[2])){
    miss_ix <- which(thresh > lpb$mx[2])
    fit_lin <- lm(log10(lambda_per_base) ~ mx, data=lpb[seq(num_points_project),])
    pred_miss <- 10^(predict(fit_lin, newdata = data.frame("mx" = thresh[miss_ix])))
    preds[miss_ix] <- pred_miss
  }
  return(preds)
}


get_thresh_with_rate <- function(rates, lpb, range = range(lpb$mx),
                                 num_points_project=10, napprox=100){
  stopifnot(all(diff(lpb$mx) <= 0))
  stopifnot(all(diff(lpb$lambda_per_base) >= 0))
  if(nrow(lpb) < num_points_project){
    warning("Fewer permutations than desired for projection. Consider running more permutations?\n")
    num_points_project <- nrow(lpb)
  }

  mxapprox <- seq(range[1], range[2], length.out=napprox)
  lpbapprox <- get_rate_with_thresh(thresh=mxapprox, lpb, num_points_project)
  preds <- approx(x=log10(lpbapprox), y=mxapprox, xout = log10(rates))$y
  #if(any(rates < lpb$lambda_per_base[2])){
  #  miss_ix <- which(rates < min(lpb$lambda_per_base))
  #  fit_lin <- lm(log10(lambda_per_base)~mx, data=lpb[seq(num_points_project),])
  #  pred_miss <- (log10(rates[miss_ix]) - fit_lin$coefficients[1])/fit_lin$coefficients[2]
  #  preds[miss_ix] <- pred_miss
  #}
  return(preds)
}
