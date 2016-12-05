collect_fret_stats <- function(temp.dir, temp.prefix){
  fl <- list.files(temp.dir, pattern = temp.prefix, full.names = TRUE)
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
    sts <- rbind(sts, R$sts)
    if(k > 1) sts.smooth <- rbind(sts.smooth, R$sts.smooth)
    if(k > 2) m1 <- rbind(m1, R$m1)
    if(k > 3){
      perm.var <-rbind(perm.var, R$perm.var)
      mperm <- rbind(mperm, R$mperm)
    }
    if(permstats) sts.perm.smooth <- rbind(sts.perm.smooth, R$sts.perm.smooth)
  }
  cat("\n")
  R$sts <- sts
  if(k > 1) R$sts.smooth <- sts.smooth
  if(k > 2) R$m1 <- m1
  if(k > 3) R$mperm <- mperm
  if(permstats) R$sts.perm.smooth <- sts.perm.smooth
  return(R)
}
