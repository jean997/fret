#'@import wavethresh


#'@title Run waveQTL
#'@description Runs the method of Shim and Stephens 2014
#'@param dat matrix or data frame of data. First column should be position.
#' second column should be window number, the rest should be data.
#'@param x Phenotype
#'@param waveQTL_loc  Path to the waveQTL program
#'@param names Names column in pheno.file
#'@return A list
#'@export
run_waveQTL <- function(dat, x,
                        waveQTL_loc="~/WaveQTL-master/bin/WaveQTL"){

  #waveQTL needs to output a lot of data -- N will be the label for these files
  N <- floor(runif(n=1, min=10, max=1e9)) #Random number label
  geno <- c("chr1.1", "A", "G", x) #fake genotype info
  cat(geno, file=paste0("geno_", N, ".txt"))

  n <- ncol(dat)
  k <- nrow(windows)

  ww <- rle(dat$win)
  starts <- c(1, cumsum(ww$lengths)+1)
  starts <- starts[-length(starts)]
  stops <- cumsum(ww$lengths)
  dat.ix <- 3:ncol(dat)
  pvals <- c()
  for(i in 1:k){
    pheno.dat <- t(dat[starts[i]:stops[i], dat.ix])
    res <- WaveQTL_preprocess(Data = pheno.dat, meanR.thresh = 0)
    f <- tempfile(tmpdir = ".")
    write.table(res$WCs, file=paste0(f, "_pheno.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
    cat(res$filtered.WCs, file=paste0(f, "_use.txt"))
    cmd <- paste0(waveQTL_loc, " -gmode 1 -g geno_", N, ".txt -p ",
            f, "_pheno.txt -u ", f, "_use.txt -o temp", N, " -f ", n, " -numPerm 100000 -fph 2")
    system(cmd)
    pval <- read.table(paste0("output/temp", N, ".fph.pval.txt"), header=TRUE)
    pvals <- c(pvals, pval[3, 1])
    unlink(paste0(f, "_pheno.txt"))
    unlink(paste0(f, "_use.txt"))
  }
  qvals <- p.adjust(pvals, method="BH")
  rates <- data.frame(t(sapply(sort(-log10(pvals)), FUN=function(xx){
    tpr_nfp(signal, discoveries=w[-log10(pvals) >= xx, , drop=FALSE])
  })))
  rates_at <- t(sapply(level, FUN=function(ll){
    tpr_nfp(signal=signal, discoveries = w[qvals<= ll, ])
  }))

  unlink(paste0("geno_", N, ".txt"))
  return(list("pvals"=pvals, "qvals"=qvals, "rates"=rates, "rates_at"=rates_at))
}
