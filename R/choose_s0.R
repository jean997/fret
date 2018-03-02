#' Select s0 using the method of Tusher, Tibshirani and Chu (2001)
#' @param beta Effect size estimates
#' @param se Estimates of sqrt(var(beta))
#'@export
choose_s0 <- function(beta, se){

  if(any(is.na(se))){
    nmiss <- sum(is.na(se))
    cat("Warning: ", nmiss, " positions have missing sd.\n")
    ix <- which(is.na(se))
    beta <- beta[-ix]
    se <- se[-ix]
  }

  salpha <- as.numeric(quantile(se, probs=seq(0.01, 0.99, by=0.01)))
  nn <- length(salpha)
  ix <- sapply(se, FUN=function(w){ sum(w > c(0, salpha))})

  fct <- function(s0, beta, se, ix){
    #cat(s0, " ")
    xx <- beta/(se + s0)
    v <- as.numeric(by(data=xx, INDICES = ix, FUN = mad, constant = 1/0.64))
    cv <- sd(v)/mean(v)
    #m <- max(abs(v-median(v))/mad(v, constant=1))
    #return(m)
    return(cv)
  }

  zz <- sapply(salpha, FUN=fct, beta=beta, se=se, ix=ix)
  return(salpha[which.min(zz)])
}
