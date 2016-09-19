#'@export
strsplit_helper <- function(list, split, field, fixed=FALSE){
  x <- unlist(lapply(list, FUN=function(x){
    unlist(strsplit(x, split, fixed=fixed))[field]}))
}

#Coppied from GWAS Tools
#' Get an object from an .RData file
#' @export
getobj <- function (Rdata){
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}

#'@export
excursions <- function(ys, z0){
  exc <- function(y, z){
    q0 <-rle( abs(y) > z )
    p0 <- length(q0$lengths)
    ivls <- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
    return(ivls)
  }
  if(length(z0)==1){
    return(exc(ys, z0))
  }else{
    stopifnot(z0[1] > 0 & z0[2] < 0)
    ys_pos <- ys; ys_pos[ys < 0] <- 0
    ivls_pos <- exc(ys_pos, z0[1])
    ys_neg <- abs(ys); ys_neg[ys > 0] <- 0
    ivls_neg <- exc(ys_neg, -z0[2])
    ivls <- rbind(ivls_pos, ivls_neg)
    return(ivls)
  }
}

#'@export
mxlist <- function(ys, z0, zmin, return.ix = FALSE){
  stopifnot(length(z0)==length(zmin))
  stopifnot(length(z0) %in% c(1, 2))

  if(length(z0)==1){
    stopifnot(z0 > 0 & zmin >= z0)
    if(all(abs(ys) < z0)) return(c())
    #Excursions at z0
    ivls <- excursions(ys, z0)
    #Max stat value inside each excurion
    mxs <- apply(ivls, MARGIN=1, FUN=function(iv){ max(abs(ys)[iv[1]:iv[2]])})
    #o <- order(abs(mxs), decreasing=FALSE)
    #mxs <- mxs[o]
    #n <- length(mxs)
    #if(sum(abs(mxs) >= zmin) > 5) ii <- which(abs(mxs) >= zmin)
    #  else ii <- 1:(min(5, n))
    #mxs <- mxs[ii]
    if(return.ix){
      ixs <- apply(ivls, MARGIN=1, FUN=function(iv){
        (iv[1]:iv[2])[which.max(abs(ys)[iv[1]:iv[2]])]
        })
      #ixs <- ixs[o][ii]
      return(data.frame("mx"=mxs, "ix"=ixs, "iv1"=ivls[,1], "iv2"=ivls[,2]))
    }
    return(mxs)
  }else{
    stopifnot(z0[1] > 0 & z0[2] < 0)
    stopifnot(zmin[1] > 0 & zmin[2] < 0)

    if(all(ys < z0[1] & ys > z0[2])) return(c())
    #Excursions at z0 - positive
    ivls <- excursions(ys, z0)
    #Max stat value inside each excurion
    mxs <- apply(ivls, MARGIN=1, FUN=function(iv){
      yy <- ys[iv[1]:iv[2]]
      return(max(abs(yy))*sign(yy[1]))
    })
    #o <- order(mxs, decreasing=TRUE)
    #mxs <- mxs[o]
    #np <- sum(mxs > 0)
    #nn <- sum(mxs < 0)
    #if(sum(mxs >=zmin[1]) > 5 & sum(mxs<=zmin[2]) > 5) ii <- which(mxs >=zmin[1]| mxs<=zmin[2])
    #  else ii <- c(1:(min(5, np)), max((np+nn-4), np+1):(np+nn))
    #mxs <- mxs[ii]
    if(return.ix){
      ixs <- apply(ivls, MARGIN=1, FUN=function(iv){
        (iv[1]:iv[2])[which.max(abs(ys)[iv[1]:iv[2]])]
      })
      #ixs <- ixs[o][ii]
      return(data.frame("mx"=mxs, "ix"=ixs, "iv1"=ivls[,1], "iv2"=ivls[,2]))
    }
    return(mxs)
  }
}

#'@export
lamtab <- function(mx, zmin, nbp, n.perm){
  stopifnot(length(zmin) %in% c(1, 2))
  mx <- sort(mx, decreasing=TRUE)
  if(length(zmin)==1){
    stopifnot(all(mx > 0))
    return(cbind(mx, (1:length(mx))/(n.perm*nbp)))
  }else{
    ct <- c()
    npos <- sum(mx > 0)
    if(npos > 0) ct <- 1:npos
    nneg <- sum(mx < 0)
    if(nneg > 0) ct <- c(ct, nneg:1)
    return(cbind(mx,  ct/(n.perm*nbp)))
  }
}
