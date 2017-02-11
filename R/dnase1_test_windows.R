#' Run a variety of tests in pre-defined windows
#'@description Calculate test statistics binning over regions
#'@param dat.file Data file
#'@param pheno.file Phenotype file
#'@param maxit Maximum iterations for Huber estimation
#'@param trait Trait name in pheno.file
#'@param names Names column in pheno.file
#' @return A list
#'@export
dnase1_test_windows <- function(dat.file, pheno.file, window.file, chr,
													maxit=50, trait="pheno", names="name",
													save.collapsed=NULL, waveqtl.min.pval =1e-8, 
                          waveQTL_loc="~/WaveQTL-master/bin/WaveQTL"){
  #normp <- function(x){
  #  return(2*pnorm(abs(x), lower.tail=FALSE))
  #}
  tp <- function(x, df){
    return(2*pt(abs(x), df=df, lower.tail=FALSE))
  }
  pois_reg <- function(y, labs, zero.val=1e-11){
    cts <- as.vector( by(y, labs, FUN=function(z){max(zero.val, sum(z))}))
    n <- as.vector( by(labs, labs, FUN=length))
    #poisson
    beta1 <- log(cts[2]/n[2])-log(cts[1]/n[1])
    mu = rep(cts[1]/n[1], sum(n))
    mu[labs==1] = cts[2]/n[2]
    phi = 1/(sum(n) - 2)* sum( (y-mu)^2/mu)
    s1 = sqrt(phi)*sqrt(sum(1/cts))
    return(c(beta1, s1, beta1/s1, tp(beta1/s1, df=sum(n)-2)))
  }
	huber_reg <- function(y, labs){
    f <- rlm(y~labs, psi=psi.huber, k=1.345, scale.est="Huber", maxit=50)
    b1 <- summary(f)$coefficients[2, 1]
    s <- summary(f)$coefficients[2, 2]
    if(is.na(s)) return(c(0, 1))
    return(c(b1, s, b1/s, tp(b1/s, df=length(y)-2)))
  }
  tt <- function(y, labs){
    b <- mean(y[labs==1]-mean(y[labs==0]))
    s <- sqrt( var(y[labs==1])/sum(labs==1) + var(y[labs==0])/sum(labs==0))
    return(c(b, s, b/s, tp(b/s, df=length(y)-2)))
  }

	#dat <- getobj(dat.file)
	dat <- read_delim(dat.file, delim=" ")
  X <- read_delim(pheno.file, delim=" ")
	datnames = names(dat)[3:ncol(dat)]
	stopifnot(all(datnames %in% X[[names]]))
  X <- X[match(datnames, X[[names]]),  ]
  x <- X[[trait]]
	win.ref <- read_delim(window.file, delim="\t", col_names=FALSE)
	win.ref <- win.ref[win.ref[,1]==chr,]

	wins <- unique(dat$win)
	stopifnot(length(wins)==nrow(win.ref))
	stopifnot(all(wins==1:nrow(win.ref)))
	dat.collapsed <- by(data=dat[, 3:ncol(dat)], INDICES=dat$win, FUN=colSums)
	dat.collapsed <- matrix(unlist(dat.collapsed), byrow=TRUE, nrow=length(wins))

  #Window start stop
  #Four collumns for each stat: beta, se, stat, p-value
  res <- apply(dat.collapsed, MARGIN=1, FUN=function(y){
    sum0 = sum(y[x==0])
    sum1 = sum(y[x==1])
    c(huber_reg(y, x), pois_reg(y, x), tt(y, x), sum0, sum1)
  })
  res <- data.frame(t(res))
	res <- cbind(win.ref, res)
  nms <- c("Chromosome", "Start", "Stop")
  nms <- c(nms, paste0("Huber_", c("Beta", "SE", "Stat", "P")))
  nms <- c(nms, paste0("Pois_", c("Beta", "SE", "Stat", "P")))
  nms <- c(nms, paste0("T_", c("Beta", "SE", "Stat", "P")))
  nms <- c(nms, "total0", "total1")
  names(res) <- nms
  cat("Running waveQTL\n")
  waveqtl_pvals <- run_waveQTL(dat, x, min.pval=min.pval, 
                        waveQTL_loc=waveQTL_loc)
  res$waveQTL_P <- waveqtl_pvals
  cat("\n")
	if(!is.null(save.collapsed)){
		dat.collapsed <- data.frame(cbind(win.ref, dat.collapsed))
		names(dat.collapsed) <- c("Chromosome", "Start", "Stop", datnames)
		save(dat.collapsed, file=save.collapsed)
	}
  return(res)
}

