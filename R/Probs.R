#' Calculate Normalized Forward Probabilities from a Forward Table Object
#'
#' Provides an easy function for calculating the normalized forward probabilities from a forward table object fwd.
#'
#' The (j,i)-th element of of the returned matrix is the probability that j is copied by i at locus fwd$l given the haplotypes
#' observed from locus 1 up through locus fwd$l.  Note: each column represents an independent HMM.
#'
#' By convention, every diagonal element is zero.
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#'
#' @param log logical; if TRUE, the forward probabilities p are returned as log(p).
#'
#' @return matrix of forward probabilities
#'
#' @seealso \code{\link{MakeForwardTable}} to generate forward table;
#'   \code{\link{Forward}} to propogate a forward table to a new locus.
#'
#' @examples
#' ForwardProbs(fwd)
#' @export ForwardProbs
ForwardProbs <- function(fwd, log = FALSE){
  if(log){
    return(sweep(log(fwd$alpha), MARGIN = 2, STATS = log(colSums(fwd$alpha)), FUN = "-"))
  }else{
    return(exp(sweep(log(fwd$alpha), MARGIN = 2, STATS = log(colSums(fwd$alpha)), FUN = "-")))
  }
}

#' Calculate Normalized Backward Probabilities from a Backward Table Object
#'
#' Provides an easy function for calculating the normalized backward probabilities from a backward table object bck.
#'
#' The (j,i)-th element of of the returned matrix is the probability of observing the haplotypes from locus bck$l+1
#' up through locus L given that j is copied by i at bck$l divided by the sum over i of those probabilities.
#' Although this matrix does not have a clean probabilistic interpretation, it captures how similar haplotypes are to each
#' other downstream of bck$l.  Note: each column represents an independent HMM.
#'
#' By convention, every diagonal element is zero.
#'
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#'
#' @param log logical; if TRUE, the backward probabilities p are returned as log(p).
#'
#' @return matrix of backward probabilities
#'
#' @seealso \code{\link{MakeBackwardTable}} to generate backward table;
#'   \code{\link{Backward}} to propogate a backward table to a new locus.
#'
#' @examples
#' BackwardProbs(bck)
#' @export BackwardProbs
BackwardProbs <- function(bck, log = FALSE){
  if(log){
    return(sweep(log(bck$beta), MARGIN = 2, STATS = log(colSums(bck$beta)), FUN = "-"))
  }else{
    return(exp(sweep(log(bck$beta), MARGIN = 2, STATS = log(colSums(bck$beta)), FUN = "-")))
  }
}


#' Calculate Posterior Marginal Probabilities from a Forward Table Object and a Backward Table Object
#'
#' Provides an easy function for calculating the posterior marginal probabilities.
#'
#' The forward and backward tables must be at the same locus in order for them to be combined to yield the posterior marginal probabilities at locus l.
#' The (i,j)-th element of of the returned matrix is the probability that j copies i at locus fwd$l=bck$l given the haplotypes observed (from locus 1 to L).
#' Note: each column represents an independent HMM.
#'
#' By convention, every diagonal element is zero.
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#'
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#'
#' @param log logical; if FALSE, the posterior marginal probabilities p are returned as log(p).
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
#' PostProbs(fwd,bck)
#' @export PostProbs
PostProbs <- function(fwd, bck, log = FALSE){
  if(fwd$l != bck$l) {
    warning("Computing dist matrix but locus position of the forward table and backward table do not match.")
  }
  if(fwd$pars.sha256 != bck$pars.sha256) {
    warning("Computing dist matrix but parameters used to calculate the forward table and backward table do not match.")
  }
  tempmat <- fwd$alpha*bck$beta
  tempmat <- sweep(log(tempmat), MARGIN = 2, STATS = log(colSums(tempmat)), FUN = "-")
  if(log){
    return(tempmat)
  }else{
    return(exp(tempmat))
  }
}


#' Calculate a Distance Matrix from a Forward Table Object and a Backward Table Object
#'
#' Provides an easy function for calculating a local distance matrix.
#'
#' The forward and backward tables must be at the same locus.
#' The (i,j)-th element of of the returned matrix is the inferred distance d_(i,j) between haplotypes j copies i at locus fwd$l=bck$l given the haplotypes observed (from locus 1 to L).
#'
#' d_(i,j) = -( log(p_(i,j)) + log(p_(j,i)) ) / 2 where p_(i,j) is the posterior marginal probability that j copies i at locus fwd$l=bck$l given the haplotypes observed (from locus 1 to L).
#'
#' By convention, d_(i,i) = 0 for all i.
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#'
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#'
#' @return matrix of distances
#'
#' @seealso
#' \code{\link{PostProb}} to calculate the posterior marginal probabilities p_(i,j);
#' \code{\link{MakeForward}} to generate Forward table;
#' \code{\link{Forward}} to propogate a Backward table to a new locus;
#' \code{\link{MakeBackwardTable}} to generate a Backward table;
#' \code{\link{Backward}} to propogate a Backward table to a new locus.
#'
#' @examples
#' DistMat(fwd,bck)
#' @export DistMat
DistMat <- function(fwd, bck){
  if(fwd$l != bck$l) {
    warning("Computing dist matrix but locus position of the forward table and backward table do not match.")
  }
  if(fwd$pars.sha256 != bck$pars.sha256) {
    warning("Computing dist matrix but parameters used to calculate the forward table and backward table do not match.")
  }
  tempmat <- fwd$alpha*bck$beta
  tempmat <- sweep(-log(tempmat), MARGIN = 2, STATS = log(colSums(tempmat)), FUN = "+")
  diag(tempmat) <- 0
  return((tempmat + t(tempmat))/2)
}