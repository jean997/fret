
#'@title Calculate lambda and fdr for each peak
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param file.list List of files with output from fret_rates_prelim
#'@param fdr.max Maximum fdr to keep data for
#'@param perm The permutation number to be treated as original data.
#'For main results, use perm=0. If you want to extract significant regions for
#'permuted data set perm to an integer in 1:n.perm
#' @return An object that can be passed to fret_thresholds
#'@export
fret_rates <- function(file.list, fdr.max=0.8, perm=0){
  #Get info from files
  R <- getobj(file.list[1])

  zmin <- R$zmin
  s <- length(zmin)
  stopifnot(s %in% c(1, 2))
  if(s==2) stopifno(zmin[1] > 0 | zmin[2] < 0 )

  if(perm==0){
    max1 <- R$m1
    max.perm <- R$mperm
  }else{
    stopifnot(perm %in% 1:R$n.perm)
    max1 <- R$mperm[R$mperm$perm==perm, ]
    max.perm <- R$mperm[R$mperm$perm!=perm,]
    if(s == 1){
      max1 <- max1[abs(max1$mx) >= R$zmin, ]
    }else{
      ix <- which(( max1$mx > 0 & max1$mx > zmin[1]) | (max1$mx < 0 & max1$mx < zmin[2]))
      max1 <- max1[ix, ]
    }
  }
  R$seg.bounds$name <- paste0(R$seg.bounds$chr, ".", 1:nrow(R$seg.bounds))
  segment.bounds <- R$seg.bounds
  n.perm <- R$n.perm
  if(perm > 0) n.prem <- R$n.perm -1

  file.list <- file.list[-1]
  for(f in file.list){
    cat(f, "\n")
    R <- getobj(f)
    stopifnot(R$zmin==zmin)
    stopifnot(R$n.perm==n.perm)
    if(perm==0){
      max1 <- rbind(max1, R$m1)
      max.perm <- rbind(max.perm, R$mperm)
    }else{
      m1 <- R$mperm[R$mperm$perm==perm, ]
      mp <- R$mperm[R$mperm$perm!=perm,]
      max.perm <- rbind(max.perm, mp)
      if(s == 1){
        m1 <- m1[abs(max1$mx) >= R$zmin, ]
      }else{
        ix <- which(( max1$mx > 0 & max1$mx > zmin[1]) | (max1$mx < 0 & max1$mx < zmin[2]))
        m1 <- m1[ix, ]
      }
      max1 <- rbind(max1, m1)
    }
    R$seg.bounds$name <- paste0(R$seg.bounds$chr, ".", 1:nrow(R$seg.bounds))
    segment.bounds <- rbind(segment.bounds, R$seg.bounds)
  }
  max1$name <- paste0(max1$chr, ".", max1$segment)
  max.perm$name <- paste0(max.perm$chr, ".", max.perm$segment)


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
              "max.perm"=max.perm, "n.perm"=n.perm, "zmin"=zmin))
}

