
#'@title Calculate lambda and fdr for each peak
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param file.list List of files with output from fret_rates_prelim
#'@param fdr.max Maximum fdr to keep data for
#' @return An object that can be passed to fret_thresholds
#'@export
fret_rates <- function(file.list, fdr.max=0.8){
  #Get info from files
  R <- getobj(file.list[1])
  max1 <- R$m1
  max.perm <- R$mperm
  zmin <- R$zmin
  R$seg.bounds$name <- paste0(R$seg.bounds$chr, ".", 1:nrow(R$seg.bounds))
  segment.bounds <- R$seg.bounds
  n.perm <- R$n.perm
  file.list <- file.list[-1]
  for(f in file.list){
    cat(f, "\n")
    R <- getobj(f)
    stopifnot(R$zmin==zmin)
    stopifnot(R$n.perm==n.perm)
    max1 <- rbind(max1, R$m1)
    max.perm <- rbind(max.perm, R$mperm)
    R$seg.bounds$name <- paste0(R$seg.bounds$chr, ".", 1:nrow(R$seg.bounds))
    segment.bounds <- rbind(segment.bounds, R$seg.bounds)
  }
  max1$name <- paste0(max1$chr, ".", max1$segment)
  max.perm$name <- paste0(max.perm$chr, ".", max.perm$segment)


  #Check inputs
  s <- length(zmin)
  stopifnot(s %in% c(1, 2))

  K <- nrow(segment.bounds)
  cat("There are ", K, " segments total.\n")

  max1 <- max1[order(max1$lambda_perbase, decreasing=FALSE), ]
  if(s==2) nbp <- cbind(segment.bounds$nbp, segment.bounds$nbp)
    else nbp <- matrix(segment.bounds$nbp, nrow=K)
  #nbp <- matrix(rep(segment.bounds$nbp, s), byrow=TRUE, nrow=s)
  mlp_ix <- grep("max_lambda", names(segment.bounds))
  max.lambda <- as.matrix(segment.bounds[,mlp_ix, drop=FALSE]*nbp)

  max1$lambda <- sapply(max1$lambda_perbase, FUN=function(r){
    sum(pmin(r*nbp, max.lambda))
  })
  max1$fdr <- max1$lambda/(1:nrow(max1))

  #Don't keep all of max1 or max.perm since it's huge
  max1 <- max1[max1$fdr <= fdr.max,]
  tot.disc <- nrow(max1) ###This is the number of discoveries at fdr.max
  lam.target <- fdr.max*tot.disc
  lam.target.pb <- lam.target/sum(segment.bounds$nbp)
  max.perm <- max.perm[max.perm$lambda_perbase <= lam.target.pb,]

  return(list("max1"=max1, "segment.bounds"=segment.bounds,
              "max.perm"=max.perm, "nbp"=nbp, "zmin"=zmin))
}


match_segments <- function(chr, pos, segment.bounds,
                           parallel=FALSE, cores=parallel::detectCores()-1){
  chroms <- unique(chr)
  #o <- order(segment.bounds$chr, segment.bounds$start, decreasing=FALSE)
  #stopifnot(all(o==1:nrow(segment.bounds)))
  ix <- rep(NA, length(chr))
  if(parallel){
    cl <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(cl))
  }
  for(c in chroms){
    cat(c, "..")
    ix_dat <- which(chr==c)
    ix_sb <- which(segment.bounds$chr==c)
    stopifnot(all(pos[ix_dat] <= max(segment.bounds$stop[ix_sb])))
    if(parallel){
      slocal <- parSapply(cl, pos[ix_dat], FUN=function(p){
        which.min(c(segment.bounds$start[ix_sb], Inf) <= p)-1
      })
    }else{
      slocal <- sapply(pos[ix_dat], FUN=function(p){
        which.min(c(segment.bounds$start[ix_sb], Inf) <= p)-1
      })
    }
    stopifnot(all(pos[ix_dat]) <= segment.bounds$stop[slocal])
    stopifnot(length(ix[ix_dat])== length(ix_sb[slocal]))
    ix[ix_dat] <- ix_sb[slocal]
  }
  return(ix)
}
