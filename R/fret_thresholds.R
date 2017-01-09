#'@title Get thresholds for a target FDR level
#'@param obj Object produced by fret_rates
#'@param target.fdr trget fdr level
#'@param stats.files Two column matrix.
#'The first column is chromosome. The second column is the corresponding stats file.
#' @return An object that can be passed to fret_thresholds
#'@export
fret_thresholds <- function(obj, target.fdr, stats.files){

  stopifnot(ncol(stats.files)==2)
  chrs <- unique(obj$max1$chr)
  stopifnot(all(chrs %in% stats.files[,1]))
  stopifnot(length(target.fdr)==1)
  if(target.fdr < min(obj$max1$fdr)) stop("The smallest possible FDR is",
                                          min(obj$max1$fdr),  ".\n")

  K <- dim(obj$max.lambda.pb)[2]
  s <- dim(obj$max.lambda.pb)[1]

  thresholds <- data.frame(matrix(nrow=K, ncol=5))
  names(thresholds) <- c("num.disc", "thresh.pos", "thresh.neg", "chrom", "file")
  #For each segment record 1) # of discoveries 2) pos threshold 3) neg threshold 4) chromosome 5) file

  tot.disc <- sum(obj$max1$fdr <= target.fdr) ###This is the number of discoveries
  #We want to draw thresholds with lambda = target.fdr*total num discoveries
  lam.target <- target.fdr*tot.disc
  tt <- get_thresh_with_rate(obj$max.perm, obj$max.lambda.pb, obj$nbp,
                             lam.target, obj$zmin)
  thresholds$thresh.pos <- tt[1,]
  if(s==1) thresholds$thresh.neg <- -tt[1,]
  else thresholds$thresh.neg <- tt[2,]

  for(j in 1:K) thresholds$num.disc[j] <- sum(obj$max1$segment[1:tot.disc]==j)
  ix <- which(thresholds$num.disc > 0)
  thresholds$chrom[ix] <- obj$max1$chr[match(ix, obj$max1$segment)]
  thresholds$file[ix] <- stats.files[match(thresholds$chrom[ix], stats.files[,1]), 2]
  discoveries <- get_discoveries(max1=obj$max1[1:tot.disc,], thresholds = thresholds)
  ret <- list("thresholds"=thresholds, "discoveries"=discoveries)
  return(ret)
}

