#'@export
collect_fret_stats <- function(temp.dir, temp.prefix, which.chunk,
                               out.file=NULL, del.temp=FALSE){
  fl_exist <- list.files(temp.dir, pattern = temp.prefix, full.names = TRUE)
  fl <- paste0(temp.dir,"/", temp.prefix, ".", which.chunk, ".RData")
  stopifnot(all(fl %in% fl_exist))
  n <- length(fl)
  stopifnot(n > 0)
  cat(n, " files to collect.\n")

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
  }
  cat("\n")
  R$sts <- sts
  if(k > 1) R$sts.smooth <- sts.smooth
  if(k > 2) R$m1 <- m1
  if(k > 3) R$mperm <- mperm
  if(is.null(out.file)) return(R)
  save(R, file=out.file)
  if(del.temp){
    for(f in fl) unlink(f)
  }
}
