#Calculate stats for one trait
stats1 <- function(Y, x, s0 = 0,  cores, stat_fun, libs = c()){
  if(cores==1){
    B <- data.frame(t(
      apply(Y, MARGIN=1, FUN=function(y){
      stat_fun(y, x, s0)
    }))) %>% dplyr::rename("beta" ="X1", "se" = "X2", "stat" = "X3")
    return(B)
  }

  #cl <- makeCluster(cores, type="FORK")
  cl <- makeCluster(cores, type="PSOCK")
  clusterExport(cl, varlist = c("stat_fun"), envir=environment())
  if(length(libs) > 0){
    txt <- paste0("library(", libs, ")")
    clusterCall(cl, fun=function(txt){sapply(txt, function(x){eval(parse(text=x))})},
                txt=txt)
  }
  B <- parApply(cl, Y, MARGIN=1, FUN=function(y, myx, s0){
    stat_fun(y, myx, s0)
  }, myx=x, s0=s0)
  stopCluster(cl)
  B <- data.frame(t(B)) %>% dplyr::rename("beta" ="X1", "se" = "X2", "stat" = "X3")
  return(B)
}

#Calculate stats for many traits (used for permutations)
stats_many <- function(Y, X, s0=0, cores, stat_fun, libs=c()){
  if(cores==1){
    B <- t(apply(Y, MARGIN=1, FUN=function(y){
      apply(X, MARGIN=2, FUN=function(x){
        stat_fun(y, x, s0)[3]
      })
    }))
    return(B)
  }

  #cl <- makeCluster(cores, type="FORK")
  cl <- makeCluster(cores, type="PSOCK")
  clusterExport(cl, varlist= c("stat_fun"), envir=environment())

  if(length(libs) > 0){
    txt <- paste0("library(", libs, ")")
    clusterCall(cl, fun=function(txt){sapply(txt, function(x){eval(parse(text=x))})},
                txt=txt)
  }

  B <- parApply(cl, X=Y, MARGIN=1, FUN=function(y, myX, s0){
    apply(myX, MARGIN=2, FUN=function(x, y, s0){
      stat_fun(y, x, s0)[3]
    }, y=y, s0=s0)
  }, myX = X, s0=s0)
  stopCluster(cl)
  B <- t(B)
  return(B)
}
