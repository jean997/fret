
#'@title Calculate lambda and fdr for each peak
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param peaks List of files with output from fret_rates_prelim
#'@param segment_info Maximum fdr to keep data for
#'@param max_fdr The permutation number to be treated as original data.
#'For main results, use perm=0. If you want to extract significant regions for
#'permuted data set perm to an integer in 1:n.perm
#' @return An object that can be passed to fret_thresholds
#'@export
fret_rates <- function(peaks, segment_info, max_fdr=0.8){
  peaks$lambda <- sapply(peaks$lambda_perbase, function(l){
    lpb <- pmin(l, segment_info$max_lambda_perbase)
    sum(lpb*segment_info$length)
  })
  o <- order(peaks$lambda, decreasing = FALSE)
  peaks <- peaks[o,]
  peaks$fdr <- peaks$lambda/seq(nrow(peaks))
  peaks <- filter(peaks, fdr <= max_fdr)
  return(peaks)
}
