


collect_from_files <- function(file.list, chrom, segment.bounds=NULL, min.length=NULL, q=0.05){
  stopifnot(length(chrom)==length(file.list))
  cc <- unique(chrom)
  if(is.null(segment.bounds)){
    cat("Collecting permutation variances and setting segment.bounds.\n")
    if(is.null(min.length)) stop("Please provide either min.length or segment.bounds.\n")
    segment.bounds <- data.frame("chrom"=c(), "start"=c(), "stop"=c())
    for(c in cc){
      fl <- file.list[chrom==c]
      R <- getobj(fl[1])
      range <- R$range
      vv <- R$perm.var
      pos <- R$stats.smooth$pos
      maxpos <- pos[length(pos)]
      for(f in fl[-1]){
        R <- getobj(f)
        stopifnot(R$range[2] > range[2] & R$range[1] > range[1])
        ix <- which.max(pos > maxpos)
        n <- nrow(R$stats.smooth)
        vv <- c(vv, R$perm.var[ix:n])
        pos <- c(pos, R$stats.smooth$pos[ix:n,1])
        maxpos <- pos[length(pos)]
      }
      sb <- find_segments(vv, pos, min.length, z0 = NULL, z = NULL, bandwidth = NULL,q = 0.05)
      sb_df <- data.frame("chrom"=rep(c, nrow(sb)), "start"=sb[,1], "stop"=sb[,2])
      segment.bounds <- rbind(segment.bounds, sb)
    }
  }

  max1 <- max.perm <- max.perm.save <-  data.frame( "mx"=c(), "chr"=c(), "pos"=c(), "iv1"=c(), "iv2"=c())
  curr.segment <- 1

  for(i in 1:length(file.list)){
    f <- file.list[i]
    cat(f, "\n")
    R <- getobj(f)
    max1 <- rbind(max1, R$max1)
    #max.perm
    mp <- rbind(max.perm.save, R$max.perm)
    segment <- apply(mp[, c("chr", "pos")], MARGIN=1, FUN=function(x){
      s <- which(segment.bounds$chr==as.character(x[1]) &
                 segment.bounds$start <= as.numeric(x[2]) &
                 segment.bounds$stop >= as.numeric(x[2]) )
      if(length(s)==0) return(NA)
      return(s)
    })
    if(i < length(chrom) & chrom[i]==chrom[i+1]){ #Don't split segments
      curr.segment <- max(segment)
      max.perm.save <- mp[segment==max(segment),]
      mp <- mp[segment < max(segment),]
      segment <- segment[segment < max(segment)]
    }else{
      max.perm.save <-  data.frame( "mx"=c(), "chr"=c(), "pos"=c(), "iv1"=c(), "iv2"=c())
    }
    if(nrow(mp)==0) next
    trim.ix <- trim_maxperm(mp$mx, segment, R$zmin)
    mp <- mp[trim.ix,]
    max.perm <- rbind(max.perm, mp)
  }

  ret <- list("max1"=max1, "max.perm"=max.perm, "zmin"=R$zmin,
            "segment.bounds"=segment.bounds)
  return(ret)
}

trim_maxperm <- function(mx, segment, zmin, nmin=10){
  s <- length(zmin)
  segs <- unique(segment)
  ix_tot <- c()
  for(seg in segs){
    ix_s <- which(segment==seg)
    my_mx <- mx[ix_s]
    if(s==1){
      ix <- which(my_mx > zmin)
      if(length(ix)==0){
        ix <- c()
      }else if(length(ix) <= nmin){
        n <- min(nmin, length(ix_s))
        ix <- order(my_mx, decreasing = TRUE)[1:n]
      }
    }else{
      #s==2
      if(sum(my_mx > 0) == 0){ 
        ixpos <- c() 
      }else if(sum(my_mx > zmin[1]) >= nmin){
        ixpos <- which(my_mx > zmin[1])
      }else{
        ixpos <- order(my_mx, decreasing=TRUE)[1:min(nmin, sum(my_mx > 0))]
      }
      if(sum(my_mx < 0) == 0){ 
        ixneg <- c() 
      }else if(sum(my_mx < zmin[2]) >= nmin){
        ixneg <- which(my_mx < zmin[2])
      }else{
        ixneg <- order(my_mx, decreasing=FALSE)[1:min(nmin, sum(my_mx < 0))]
      }
      ix <- c(ixpos, ixneg)
    }
    ix_tot <- c(ix_tot, ix_s[ix])
  }
  return(ix_tot)
}
