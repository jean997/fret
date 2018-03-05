
#'@export
pre_smooth <- function(file_name, bandwidth, out_file_name, maxzero,
                       windows_bed_file=NULL, chr=NULL,
                       smoother=c("ksmooth_0", "ksmooth"), chunksize=1e5,
                       out_by = floor(bandwidth/2)){

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
  if(!is.null(windows_bed_file)){
    if(is.null(chr)) stop("If you provide a bed file, please also provide the chromosome.\n")
    if(is.null(out_by)) stop("If you provide a bed file, please also provide out_by.\n")
    windows <- read_tsv(windows_bed_file, col_names=c("chrom", "start", "stop")) %>%
                filter(chrom==chr)
    new_pos_full <- apply(windows, 1, function(x){seq(x[2], x[3], by=out_by)}) %>% unlist()
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
  options(scipen=1000)
  for(i in 1:nchunks){
    goto(df_laf, max(1, (i-1)*chunksize-bandwidth + 1))
    dat <- next_block(df_laf, nrows=chunksize + 2*bandwidth)
    read_start <- max(1, (i-1)*chunksize-bandwidth + 1)

    pos <- dat$pos
    if(i==1) first_pos <- pos[1]
    keep_start <- first_pos
    if(i > 1) keep_start <- first_pos + 1
    keep_stop <- pos[chunk_stop[i] - read_start + 1]
    if(!is.null(windows_bed_file)){
      pos_new <- new_pos_full[ new_pos_full >= keep_start & new_pos_full <= keep_stop]
    }else if(!is.null(out_by)){
      pos_new <- seq(first_pos, max(pos), by=out_by)
      pos_new <- pos_new[pos_new >=min(pos) & pos_new <= max(pos)]
    }else{
      pos_new <- pos
    }
    dat_sm <- apply(dat[,-1], 2, function(y){smooth_func(pos, y, pos_new, bandwidth)})

    dat_sm <- cbind(pos_new, dat_sm)  %>% as.data.frame() %>%
      rename("pos" = "pos_new") %>%
      filter(pos >= keep_start & pos <= keep_stop)
    first_pos <- dat_sm$pos[nrow(dat_sm)]
    nzero <- select(dat_sm, -pos) %>% apply(., 1, function(x){sum(x==0)})
    dat_sm <- dat_sm %>% filter(nzero <= maxzero)
    write.table(dat_sm, file=out_file_name, sep=" ", row.names=FALSE, quote=FALSE, col.names=(i==1), append = (i > 1))
    #if(17473921 %in% pos_new) break
  }
  close(df_laf)

}
