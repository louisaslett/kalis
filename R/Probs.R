#' Calculate Forward Probabilities from a Forward Table Object
#'
#' Provides an easy function for calculating the log forward probabilities from a forward table object fwd.
#'
#' The (i,j)-th element of of the returned matrix is the probability that j copies i at locus fwd$l and the haplotypes
#' observed from locus 1 up through locus fwd$l.
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#'
#' @param log logical; if TRUE (default), the forward probabilities p are returned as log(p).
#'
#' @return matrix of forward probabilities
#'
#' @seealso \code{\link{MakeForwardTable}} to generate forward table;
#'   \code{\link{Forward}} to propogate a forward table to a new locus.
#'
#' @examples
#' fwd <- MakeForwardTable()
#' Forward(fwd, 100, Pi, mu, rho)
#' ForwardProbs(fwd)
#' @export ForwardProbs
ForwardProbs <- function(fwd, log=TRUE){
  if(log==TRUE){
    return(sweep(log(fwd$alpha),2,fwd$alpha.f2))
  }else{
    return(exp(sweep(log(fwd$alpha),2,fwd$alpha.f2)))
  }
}

#' Calculate Backward Probabilities from a Backward Table Object
#'
#' Provides an easy function for calculating the log backward probabilities from a backward table object bck.
#'
#' The (i,j)-th element of of the returned matrix is the probability that j copies i at locus fwd$l and the haplotypes
#' observed from locus fwd$l+1 up through locus L.
#'
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#'
#' @param log logical; if TRUE (default), the backward probabilities p are returned as log(p).
#'
#' @return matrix of backward probabilities
#'
#' @seealso \code{\link{MakeBackwardTable}} to generate backward table;
#'   \code{\link{Backward}} to propogate a backward table to a new locus.
#'
#' @examples
#' bck <- MakeBackwardTable()
#' Backward(bck, 100, Pi, mu, rho)
#' BackwardProbs(bck)
#' @export ForwardProbs
BackwardProbs <- function(bck, log=TRUE){
  if(log==TRUE){
    return(sweep(log(bck$beta),2,bck$beta.g2))
  }else{
    return(exp(sweep(log(bck$beta),2,bck$beta.g2)))
  }
}


#' Calculate Posterior Marginal Probabilities from a Forward Table Object and a Backward Table Object
#'
#' Provides an easy function for calculating the log posterior marginal probabilities.
#'
#' Note that the forward and backward tables must be at the same locus in order for them to be combined to yield the posterior marginal probabilities.
#' The (i,j)-th element of of the returned matrix is the probability that j copies i at locus fwd$l=bck$l given the haplotypes observed (from locus 1 to L).
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#'
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#'
#' @param log logical; if TRUE (default), the posterior marginal probabilities p are returned as log(p).
#'
#' @return matrix of posterior marginal probabilities
#'
#' @seealso
#' \code{\link{MakeForward}} to generate Forward table;
#'   \code{\link{Forward}} to propogate a Backward table to a new locus;
#' \code{\link{MakeBackwardTable}} to generate a Backward table;
#'   \code{\link{Backward}} to propogate a Backward table to a new locus.
#'
#' @examples
#' fwd <- MakeForwardTable()
#' Forward(fwd, 100, Pi, mu, rho)
#' bck <- MakeBackwardTable()
#' Backward(bck, 100, Pi, mu, rho)
#' BackwardProbs(bck)
#' PostProbs(fwd,bck)
#' @export PostProbs
PostProbs <- function(fwd, bck, log=TRUE){
  if(log==TRUE){
    return(sweep(log(fwd$alpha*bck$beta),2,fwd$alpha.f2+bck$beta.g2))
  }else{
    return(exp(sweep(log(fwd$alpha*bck$beta),2,fwd$alpha.f2+bck$beta.g2)))
  }
}


