collect_fret_stats <- function(temp.dir, temp.prefix){
  fl <- list.files(temp.dir, pattern = temp.prefix)
  n <- length(fl)
  stopifnot(n > 0)
  cat(n, " files to collect.\n")

  chunks <- as.numeric(strsplit_helper(fl, ".", 2, fixed=TRUE))
  o <- order(chunks, decreasing = FALSE)
  fl <- fl[o]

  R <- getobj(fl[1])
  sts <- R$sts
  k <- 1
  if("sts.smooth" %in% names(R)){
    sts.smooth <- R$sts.smooth
    k <- k+1
  }
  if("sts.perm.smooth" %in% nanmes(R)){
    sts.perm.smooth <- R$sts.perm.smooth
    k <- k + 1
  }
  for(f in fl[-1]){
    cat(fl, "..")
    R <- getobj{f}
    sts <- rbind(sts, R$sts)
    if(k > 1) sts.smooth <- rbind(sts.perm, R$sts.perm)
    if(k > 2) sts.perm.smooth <- rbind(sts.perm.smooth, R$sts.perm.smooth)
  }
  cat("\n")
  R$sts <- sts
  if(k > 1){
    R$sts.smooth <- sts.smooth
    R$max1 <- mxlist(sts.smooth$ys, R$z0, R$zmin, return.ix = TRUE)
    if(length(R$zmin)==1) R$max1 <- R$max1[R$max1$mx > R$zmin,] #Unsigned
      else max1 <- R$max1[R$max1$mx > R$zmin[1] | R$max1$mx < R$zmin[2],] #signed
    #Convert index to position
    R$max1$pos <- sts$pos[max1$ix]
    R$max1$iv1 <- sts$pos[max1$iv1]
    R$max1$iv2 <- sts$pos[max1$iv2]
    R$max1$chr <- R$chrom
    R$max1 <- R$max1[, c("mx", "chr", "pos", "iv1", "iv2")]
    if(!is.null(R$range)) R$max1 <- R$max1[R$max1$mx >= R$range[1] & R$max1$mx <= R$range[2],]
  }
  if(k > 2){
    R$sts.perm.smooth <- sts.perm.smooth
    R$perm.var <- apply(stats.perm.smooth, MARGIN=1, FUN=var)
    max.perm <- apply(stats.perm.smooth, MARGIN=2, FUN=function(yys){
      mxlist(yys, R$z0, R$zmin, return.ix = TRUE)
    })
    R$max.perm <- do.call(rbind, max.perm)
  }
  return(R)
}
