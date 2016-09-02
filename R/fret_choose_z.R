
#' Find thresholds for a range of lambda values. Calculate FDR.
#'@description Find thresholds for a range of lambda values. Calculate FDR.
#'@param max1.list List of peak heights for each segment.
#'@param perm.maxes.list a list of matrices giving peak heights for permutations. Each matrix has two columns. First column is the threshold level,
#'second column is lambda per bp (num of peaks less than threshold divided by number of permutations and number of base-pairs).
#'@param nbp Number of base-pairs per segment
#'@param zmin Lower bound on significance thresholds.
#'@param nlam Number of lambda values to consider
#' @return A list with items z, Robs, and fdr.
#'@export
fret_choose_z <- function(max1.list, perm.maxes.list, nbp, zmin, nlam, seg.names=NULL,
                          log.lambda.min=NULL, log.lambda.max=NULL){
  stopifnot(length(perm.maxes.list)==length(nbp) & length(perm.maxes.list)==length(max1.list))
  n.seg <- length(max1.list)

  s <- length(zmin)

  if(is.null(log.lambda.min)){
    log.lambda.min <- Inf
    log.lambda.max <- -Inf
    for(i in 1:n.seg){
      m <- perm.maxes.list[[i]]
      log.lambda.min <- min(log.lambda.min, log10(s*m[,2][m[,2] > 0]))
      log.lambda.max <- max(log.lambda.max, log10(s*m[,2][m[,2] > 0]))
    }
    lams <- seq(log.lambda.min, log.lambda.max, length.out=nlam)
  }else{
    lams <- seq(log.lambda.min, log.lambda.max, length.out=nlam)
  }

  if(is.null(seg.names)) seg.names <- paste0("segment", 1:n.seg)

  #lambda per bp, lambda target, segments, num per segment
  Robs <- z <- matrix(nrow=nlam, ncol=n.seg +2 )
  z[,1] <- Robs[, 1] <- lams
  z[,2] <- Robs[,2] <- lams + log10(sum(nbp))

  if(length(zmin)==1){
    for(i in 1:n.seg){
      m <- perm.maxes.list[[i]]
      #No permutation peaks ever reached zmin --> zmin is the threshold for all lambda
      if(nrow(m)==1 & m[1, 2]==0){
        z[, i+2] = zmin
        Robs[, i+2] = sum(max1.list[[i]] >= zmin)
      }else{
        z[, i+2] <- approx(y=m[,1], x=log10(m[,2]),
                            xout=lams, yright=zmin, yleft=Inf)$y
        z[, i+2] <- pmax(zmin, z[, i+2])
        Robs[, i+2] <- sapply(z[, i+2], FUN=function(zz){sum(max1.list[[i]] > zz)})
      }
    }
    z <- data.frame(z)
    Robs <- data.frame(Robs)
    names(z)<- names(Robs) <- c("log_lambda_perbase", "log_lambda",  seg.names)
    fdr <- (10^(Robs$log_lambda))/rowSums(Robs[, c(-1, -2),  drop=FALSE])
    ret <- list("Robs"=Robs, "z"=z, "fdr"=fdr, "nbp"=nbp)
    return(ret)
  }else{
    zneg <- z
    for(i in 1:n.seg){
      m <- perm.maxes.list[[i]]
      #No permutation peaks ever reached zmin --> zmin is the threshold for all lambda
      if(nrow(m)==1 & m[1, 2]==0 | nrow(m)==0){
        z[, i+2] = zmin
        Robs[, i+2] = sum(max1.list[[i]] >= zmin[1] | max1.list[[i]] <= zmin[2])
      }else{
        npos <- sum(m[,1] > 0)
        ntot <- nrow(m)
        z[, i+2] <- approx(y=m[1:npos,1], x=log10(m[1:npos,2]),
                           xout=lams-log10(2), yright=zmin[1], yleft=Inf)$y
        z[, i+2] <- pmax(zmin[1], z[, i+2])
        zneg[, i+2] <- approx(y=m[(npos+1):ntot,1], x=log10(m[(npos+1):ntot,2]),
                           xout=lams-log10(2), yright=zmin[2], yleft=-Inf)$y
        zneg[, i+2] <- pmin(zmin[2], zneg[, i+2])
        Robs[, i+2] <- sapply(1:nlam, FUN=function(j){
          sum(max1.list[[i]] >= z[j, i+2] | max1.list[[i]] <= zneg[j, i+2])})
      }
    }
    z <- data.frame(z)
    zneg <- data.frame(zneg)
    Robs <- data.frame(Robs)
    names(z)<- names(Robs) <- names(zneg) <- c("log_lambda_perbase", "log_lambda",  seg.names)
    fdr <- (10^(Robs$log_lambda))/rowSums(Robs[, c(-1, -2),  drop=FALSE])
    ret <- list("Robs"=Robs, "z"=z, "zneg" = zneg, "fdr"=fdr, "nbp"=nbp)
    return(ret)
  }
}
