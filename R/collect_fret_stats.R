#' @title Collect fret_stats temporary files
#' @description Combine temporary files from fret_stats into one file.
#' Also determine interval boundaries.
#' @param file_list List of files produced by fret_stats
#' @param seg_type One of "by_file", "by_chromosome", "find", "provided". See description.
#' @param chromosome Vector the same length as file_list giving the chromosome for each file.
#' This is only required if seg_type is not equal to by_file.
#' @param segment_bounds Optional data frame with columns start, stop and chrom.
#' @param min_segment_length Used if seg_type="find". Minimum segment length.
#' Defaults to 50*bandwidth if not provided.
#'@export
stats_to_rates <- function(file_list, seg_type=c("by_file", "by_chromosome", "find", "provided"),
                           chromosome, segment_bounds,
                           min_segment_length=NULL){

  #Options
  seg_type <- match.arg(seg_type)
  if(seg_type!="by_file" & missing(chromosome)) stop("chromosome must be provided unless seg_type=\"by_file\".\n")
  R <- readRDS(file_list[1])
  if(R$n_perm==0){
    cat("There are no permutations. Statistics will be collected and returned")
    if(missing(chromosome)) chromosome <- rep("missing", length(file_list))
    chrs <- unique(chromosome)
    stats <- lapply(chrs, function(c){
      dat <- lapply(file_list[chromosome==c], function(file){
        f <- readRDS(file)
        f$stats$label <- f$label
        return(f$stats)
      })
      sts <- do.call(rbind, dat)
      sts$chromosome <- c
      return(sts)
    })
    stats <- do.call(rbind, stats)
    return(list("stats"=stats))
  }


  if(!missing(segment_bounds) & !seg_type=="find") stop("segment_bounds is only required for seg_type = \"find\".\n")
  if(seg_type=="find"){
    if(is.null(min_segment_length)){
      min_segment_length <- 50*R$bandwidth
    }
    cat("Segments will be determined from the data. Minimum segment length: ", min_segment_length, "\n")
  }
  if(seg_type=="provided"){
    if(missing(segment_bounds)) stop("If seg_type=\"provided\", segment_bounds must be provided.\n")
    stopifnot(all(c("start", "stop", "chrom") %in% names(segment_bounds)))
  }

  n <- length(file_list)
  if(!missing(chromosome)) stopifnot(length(chromosome)==n)
  cat("Statistics in ", n, " files will be collected.\n")

  if(seg_type=="by_file"){
    rate_info <- lapply(seq(file_list), function(i){
                    rts <- fret_rates_prelim(file_list[i])
                    rts$segment_info$segment <- i
                    rts$peaks <- rts$peaks %>% mutate(segment = i)
                    rts$stats$file <- i
                    rts
                  })
      peaks <- do.call(rbind, lapply(rate_info, function(x){x$peaks}))
      segment_info <- do.call(rbind, lapply(rate_info, function(x){x$segment_info}))
      stats <- do.call(rbind, lapply(rate_info, function(x){x$stats}))
      peaks <- fret_rates(peaks, segment_info)
      return(list("peaks"=peaks, "segment_info" = segment_info, "stats"=stats))
  }

  if(seg_type=="by_chromosome"){
    chrs <- unique(chromosome)
    rate_info <- lapply(chrs, function(c){
      rts <- fret_rates_prelim(file_list[chromosome==c])
      rts$segment_info$segment <- c
      rts$peaks <-  rts$peaks %>%  mutate(segment = c)
      rts$stats$chromosome <- c
      rts
    })
    peaks <- do.call(rbind, lapply(rate_info, function(x){x$peaks}))
    segment_info <- do.call(rbind, lapply(rate_info, function(x){x$segment_info}))
    stats <- do.call(rbind, lapply(rate_info, function(x){x$stats}))
    peaks <- fret_rates(peaks, segment_info)
    return(list("peaks"=peaks, "segment_info" = segment_info, "stats"=stats))
  }
  chrs <- unique(chromosome)
  if(seg_type=="find"){
    sbs <- lapply(chrs, function(c){
      fl <- file_list[chromosome==c]
      dat <- lapply(fl, function(file){
        f <- readRDS(file)
        f$stats[, c("pos", "perm_var")]
      })
      dat <- do.call(rbind, dat)
      sb <- find_segments(dat[,2], dat[,1], min_segment_length)
      sb$chrom <- c
      sb
    })
    segment_bounds <- do.call(rbind, sbs)
  }
  rate_info <- lapply(chrs, function(c){
    sb <- dplyr::filter(segment_bounds, chrom==c)
    sb <- sb[order(sb$start),]
    rts <- fret_rates_prelim(file_list[chromosome==c], segment_bounds=sb)
    rts$segment_info$chrom <- c
    rts$segment_info$segment <- paste0(c, "-", rts$segment_info$segment)
    rts$peaks <- rts$peaks %>% mutate(segment = paste0(c, "-", segment))
    rts$stats$chromosome <- c
    rts
  })
  peaks <- do.call(rbind, lapply(rate_info, function(x){x$peaks}))
  segment_info <- do.call(rbind, lapply(rate_info, function(x){x$segment_info}))
  stats <- do.call(rbind, lapply(rate_info, function(x){x$stats}))
  peaks <- fret_rates(peaks, segment_info)
  return(list("peaks"=peaks, "segment_info" = segment_info, "stats"=stats))
}

