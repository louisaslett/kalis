##### Gold Master #####

goldmaster.blobby <- function(alpha,beta,recipient_hap, unit.dist = 1, thresh = 0.2){
  f <- function(x,c){ifelse(x<c,0,pmin(x,1))}
  probs <- alpha * beta
  d <- -log(probs/sum(probs))
  d[recipient_hap] <- 0
  d[!is.finite(d)] <- 744.4400719213812180897
  sigma <- order(d)
  psi <- f(c(diff(d[sigma])/unit.dist,0),thresh)/seq_len(length(alpha))
  psi[sigma] <- rev(cumsum(rev(psi)))
  psi
}


goldmaster.blobby.dedip <- function(alpha,beta,left_recipient_hap, unit.dist = 1, thresh = 0.2){
  if(ncol(alpha) != 2){stop("alpha must be a matrix with 2 columns")}
  if(ncol(beta) != 2){stop("beta must be a matrix with 2 columns")}

  if(left_recipient_hap <= 0 || left_recipient_hap > nrow(alpha)-1 || as.integer(left_recipient_hap)!=left_recipient_hap){
    stop("left_recipient_hap must be an integer in [1,nrow(alpha)-1]")
  }
  if(c < 0 || c > 1){stop("c must be in [0,1]")}

  v <- goldmaster.blobby(alpha[,1],beta[,1],recipient_hap=left_recipient_hap,unit.dist,thresh)
  v <- v + goldmaster.blobby(alpha[,2],beta[,2],recipient_hap=left_recipient_hap+1L,unit.dist,thresh)
  v[seq.int(1,length(v),2)] + v[seq.int(2,length(v),2)]
}

goldmaster.blobby.full <- function(alpha,beta,left_recipient_hap, unit.dist, thresh){
  if(ncol(alpha) != ncol(beta)){stop("alpha and beta must have the same number of columns")}
  if(nrow(alpha) != nrow(beta)){stop("alpha and beta must have the same number of rows")}

  if( nrow(alpha)%%2 || ncol(alpha)%%2 ){stop("alpha must be a matrix with an even number of rows and columns")}

  if(left_recipient_hap <= 0 || left_recipient_hap > nrow(alpha)-1 || as.integer(left_recipient_hap)!=left_recipient_hap){
    stop("left_recipient_hap must be an integer in [1,nrow(alpha)-1]")
  }

  if(thresh < 0 || thresh > 1){stop("thresh must be in [0,1]")}

  n.samps <- ncol(alpha)/2

  M <- matrix(0,n.samps,n.samps)

  for(i in 1:n.samps){
    v <- goldmaster.blobby(alpha[,2*i-1],beta[,2*i-1],recipient_hap=left_recipient_hap+2*i-2L,unit.dist,thresh)
    v <- v + goldmaster.blobby(alpha[,2*i],beta[,2*i],recipient_hap=left_recipient_hap+2*i-1L,unit.dist,thresh)
    M[,i] <- v[seq.int(1,length(v),2)] + v[seq.int(2,length(v),2)]
  }
  M
}


CladeMat.GM <- function(fwd,bck,unit.dist,thresh = 0.2){
  M <- goldmaster.blobby.full(fwd$alpha,bck$beta,left_recipient_hap = bck$from_recipient,unit.dist,thresh)
  M
}