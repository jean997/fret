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


mxlist <- function(ys, z0, zmin){
  stopifnot(length(z0)==length(zmin))
  stopifnot(length(z0) %in% c(1, 2))

  if(length(z0)==1){
    stopifnot(z0 > 0 & zmin >= z0)
    if(all(abs(ys) < z0)) return(c())
    q0 <-rle( abs(ys) > z0 )
    p0 <- length(q0$lengths)
    #Excursions at z0
    ivls <- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
    #Max stat value inside each excurion
    mxs <- apply(ivls, MARGIN=1, FUN=function(iv){ max(abs(ys)[iv[1]:iv[2]])})
    mxs <- mxs[mxs >= zmin]
    return(mxs)
  }else{
    sgn <- z0_l <- sign(ys)
    z0_l <- replace(z0_l, sgn==1, z0[1])
    z0_l <- replace(z0_l, sgn==-1, z0[2])
    if(all(abs(ys) < abs(z0_l))) return(c())
    q0 <-rle( abs(ys) > abs(z0_l) )
    p0 <- length(q0$lengths)
    #Excursions at z0
    ivls <- cbind(c(1, cumsum(q0$lengths)[-p0]+1)[q0$values], (cumsum(q0$lengths))[q0$values])
    #Max stat value inside each excurion
    mxs <- apply(ivls, MARGIN=1, FUN=function(iv){
      yy <- ys[iv[1]:iv[2]]
      return(max(abs(yy))*sgn(yy[1]))
    })
    mxs <- mxs[mxs >= zmin[1] | mxs <= zmin[2]]
    return(mxs)
  }
}

