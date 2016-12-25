#' @title Collect fret_stats temporary files
#' @description Combine temporary files from fret_stats into one file.
#' Also determine interval boundaries.
#' @param temp.dir Directory containing temporary files
#' @param temp.prefix Prefix of temporary files (e.g. If your files are named
#' test_chr1.N.RData where N is the chunk number then temp.prefix is test_chr1)
#' @param which.chunk Vector of integers. Which chunks to collect
#' @param out.file Name an output file
#' @param del.temp Delete the temporary files when done?
#' @param min.interval.width Minimum interval length when automatically determining interval endpoints
#'@export
collect_fret_stats <- function(temp.dir, temp.prefix, which.chunk,
                               out.file=NULL, del.temp=FALSE, min.interval.width=NULL){
  fl_exist <- list.files(temp.dir, pattern = temp.prefix, full.names = TRUE)
  fl <- paste0(temp.dir,"/", temp.prefix, ".", which.chunk, ".RData")
  stopifnot(all(fl %in% fl_exist))
  n <- length(fl)
  stopifnot(n > 0)
  cat(n, " files to collect.\n")

  R <- getobj(fl[1])
  sts <- R$sts
  chrom <- R$chrom
  zmin <- R$zmin
  k <- 1
  if("sts.smooth" %in% names(R)){
    sts.smooth <- R$sts.smooth
    k <- k+1
  }
  if("m1" %in% names(R)){
    m1 <- R$m1
    k <- k+1
  }
  if("mperm" %in% names(R)){
    perm.var <- R$perm.var
    mperm <- R$mperm

    k <- k+1
  }
  for(f in fl[-1]){
    cat(f, "..")
    R <- getobj(f)
    stopifnot(R$chrom == chrom)
    stopifnot(R$zmin == zmin)
    sts <- rbind(sts, R$sts)
    if(k > 1) sts.smooth <- rbind(sts.smooth, R$sts.smooth)
    if(k > 2) m1 <- rbind(m1, R$m1)
    if(k > 3){
      perm.var <-rbind(perm.var, R$perm.var)
      mperm <- rbind(mperm, R$mperm)
    }
  }
  if(k > 2) m1$chr <- chrom
  if(k > 3) mperm$chr <- chrom
  cat("\n")
  R$sts <- sts
  if(k > 1) R$sts.smooth <- sts.smooth
  if(k > 2) R$m1 <- m1
  if(k > 3){
    R$mperm <- mperm
    R$perm.var <- perm.var
  }
  if(is.null(min.interval.width)) min.interval.width=50*R$bandwidth
  sb <- find_segments(vv=R$perm.var$var, pos=R$perm.var$pos, min.length=min.interval.width)
  R$seg.bounds <- data.frame("chr"=rep(R$chrom, nrow(sb)), "start"=sb[,1], "stop"=sb[,2])

  if(is.null(out.file)) return(R)
  save(R, file=out.file)
  if(del.temp){
    for(f in fl) unlink(f)
  }
}
