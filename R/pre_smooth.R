
#'@export
pre_smooth <- function(file_name, bandwidth, out_file_name, maxzero,
                       smoother=c("ksmooth_0", "ksmooth"), chunksize=1e5, out_by = floor(bandwidth/2)){

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

  #Open pheno file
  dm <- detect_dm_csv(filename=file_name, header=TRUE, sep=" ")
  dm$columns$type <- rep(c("integer", "double"), c(1, nrow(dm$columns)-1))
  df_laf <- laf_open(dm)

  #Determine number of chunks
  nl <- determine_nlines(file_name)-1
  nchunks <- ceiling(nl/chunksize)
  chunk_start <- (seq(nchunks) -1)*chunksize + 1
  if(nchunks == 1) chunk_stop <- nl
  	else chunk_stop <- c(chunk_start[-1] -1,  nl)
  options(scipen=Inf)
  for(i in 1:nchunks){
    goto(df_laf, max(1, (i-1)*chunksize-bandwidth + 1))
    dat <- next_block(df_laf, nrows=chunksize + 2*bandwidth)
    read_start <- max(1, (i-1)*chunksize-bandwidth + 1)
    keep_start <- chunk_start[i] -read_start + 1
    keep_stop <- chunk_stop[i] - read_start + 1
    pos <- dat$pos
    if(!is.null(out_by)) pos_new <- seq(min(pos), max(pos), by=out_by)
      else pos_new <- pos
    dat_sm <- apply(dat[,-1], 2, function(y){smooth_func(pos, y, pos_new, bandwidth)})
    nzero <- apply(dat_sm, 1, function(x){sum(x==0)})
    dat_sm <- cbind(pos_new, dat_sm)  %>% as.data.frame() %>% filter(nzero <= maxzero)
    write.table(dat_sm, file=out_file_name, sep=" ", row.names=FALSE, quote=FALSE, col.names=(i==1), append = (i > 1))
  }
  close(df_laf)

}
