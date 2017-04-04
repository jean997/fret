#'@export
qqthin <- function(pvals, thin=c(0.25, 100), shade=TRUE, truncate=Inf, group=NULL){
  n <- length(pvals)
  pvals <- sort(pvals)
  v = qunif(p=seq(0, 1, length.out = n+1)[2:(n+1)])
  v = -log10(v)

  if(!is.null(thin)){
    q = quantile(v, probs=thin[1])
    q_idx = min(which(v <= q))
    thin_idx = unique(ceiling(seq(n-q_idx, n, length.out=thin[2])))
    thin_idx = c(1:(n-q_idx))
    v=v[thin_idx]
    pvals = -log10(pvals[thin_idx])
    if(!is.null(group)) group <- group[thin_idx]
  }else{
    thin_idx <- 1:n
    pvals <- -log10(pvals)
  }
  #Shading
  if(shade){
    c975 <- sapply(thin_idx, FUN=function(i){
      qbeta(0.975, i, n - i + 1)
    })
    c025 <- sapply(thin_idx, FUN=function(i){
      qbeta(0.025, i, n - i + 1)
    })
    df.shade = data.frame("x"=c(v, rev(v)), "y"=c(-log10(c025), rev(-log10(c975))))
  }

  if(is.finite(truncate) & any(pvals > truncate)){
    shape <- rep(1, length(pvals))
    shape[pvals > truncate] <- 2
    pvals <- pmin(truncate, pvals)
  }else{
    shape=rep(1, length(pvals))
  }


  df <- data.frame("pval"=pvals, "v"=v, "shape"=shape)
  if(!is.null(group)) df$group <- group
  h <- ggplot(df)
  if(shade) h <-  h + geom_polygon(data=df.shade, aes(x=x, y=y), alpha=0.3, fill="black")
  if(!is.null(group)) h <- h +  geom_point(aes(x=v, y=pval, color=group, shape=as.factor(shape)))
    else h <- h +  geom_point(aes(x=v, y=pval, shape=as.factor(shape)))
  h <- h +
    geom_abline(slope=1, intercept=0) +
    xlab(expression(Expected~~-log[10](italic( p )-value)))+
    ylab(expression(Observed~~-log[10](italic( p )-value)))+
    theme_bw(12) + theme(panel.grid=element_blank(), legend.position="none")
  return(h)
}
