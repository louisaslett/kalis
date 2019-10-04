#' Posterior marginal probabilities
#'
#' Calculated the posterior marginal probabilities at a given locus using
#' forward and backward tables propagated to that position.
#'
#' The forward and backward tables must be at the same locus in order for them
#' to be combined to yield the posterior marginal probabilities at locus
#' \eqn{l}.
#' The \eqn{(j,i)}-th element of of the returned matrix is the probability that
#' \eqn{j} is copied by \eqn{i} at the current locus, \eqn{l}, of the two
#' tables, given the haplotypes observed (over the whole sequence).
#'
#' Note that each column represents an independent HMM.
#'
#' By convention, every diagonal element is zero.
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#' @param log logical; defaults to \code{FALSE}.
#'   If \code{TRUE}, the logarithm of the posterior marginal probabilities are
#'   returned.
#'
#' @return
#'   Matrix of posterior marginal probabilities.
#'   The \eqn{(j,i)}-th element of of the returned matrix is the probability
#'   that \eqn{j} is copied by \eqn{i} at the current locus, \eqn{l}, of the two
#'   tables, given the haplotypes observed (over the whole sequence).
#'
#' @seealso
#'   \code{\link{DistMat}} to generate calculate \eqn{d_{ji}}{d_(j,i)} distances
#'   directly;
#'   \code{\link{Forward}} to propogate a Forward table to a new locus;
#'   \code{\link{Backward}} to propogate a Backward table to a new locus.
#'
#' @examples
#' \dontrun{
#' # To get the posterior probabilities at, say, locus 10 ...
#' Forward(fwd, pars, 10, nthreads = 8)
#' Backward(bck, pars, 10, nthreads = 8)
#' p <- PostProbs(fwd, bck)
#' }
#'
#' @export
PostProbs <- function(fwd, bck, log = FALSE) {
  if(fwd$l != bck$l) {
    stop("locus position of the forward table and backward table do not match.")
  }
  if(fwd$pars.sha256 != bck$pars.sha256) {
    stop("parameters by the forward table and backward table do not match.")
  }
  tempmat <- fwd$alpha*bck$beta
  tempmat <- sweep(log(tempmat), MARGIN = 2, STATS = log(colSums(tempmat)), FUN = "-")
  if(log) {
    return(tempmat)
  } else {
    return(exp(tempmat))
  }
}

#' Distance matrix
#'
#' Utility for calculating posterior marginal probabilities and distance matrices at, in between, or excluding variants.
#'
#' This computes a local probability or distance matrix based on the forward and backward
#' tables at a certain locus.  The default usage is provide forward and backward tables at the same locus \eqn{l}.
#' so that the \eqn{(j,i)}-th element of of the returned matrix is the inferred distance
#' \eqn{d_{ji}}{d_(j,i)} between haplotypes \eqn{j} and \eqn{i} at the current
#' locus, \eqn{l}, of the two tables given the haplotypes observed (over the
#' whole sequence).
#'
#' In particular,
#'
#' \deqn{d_{ji} = -\frac{log(p_{ji}) + log(p_{ij})}{2}}{d_(j,i) = -( log(p_(j,i)) + log(p_(i,j)) ) / 2}
#'
#' where \eqn{p_{ji}}{p_(j,i)} is the posterior marginal probability that
#' \eqn{j} is coped by \eqn{i} at the current locus of the two tables, \eqn{l},
#' given the haplotypes observed (over the whole sequence).
#'
#' By convention, \eqn{d_{ii} = 0}{d_(i,i) = 0} for all \eqn{i}.
#'
#' A matrix of the raw posterior marginal probabilities \eqn{p_{ji}}{p_(j,i)} may be obtained by setting \code{neg.log = FALSE}, \code{standardize = FALSE}, and \code{symmetrize = FALSE}.
#'
#' This function also allows users to calculate distance matrices in between variants and also to calculate matrices that exclude a set of consecutive variants by passing a
#' backward table in beta.theta space.  If in beta.theta space, \code{bck$l} may be greater than but not equal to \code{fwd$l}.  When passed a backward table in beta.theta space, the transition probabilities
#' \code{rho} between \code{fwd$l} and \code{bck$l} will be summed together to form a total, \code{rho_total}.  If \code{bias} is provided as named element in the list \code{beta.theta.opts}, then the tables will be propagated to meet in between with the given \code{bias}.  More explicitly, the forward table will be propagated \code{rho_total * bias} and
#' the backward table will be propagated a distance of \code{rho_total * (1-bias)}.  The distance to propagate the forward and backward tables respectively can be given by passing \code{rho1step.fwd} and \code{rho1step.bck} as arguments in \code{beta.theta.opts}, which override \code{bias}.
#'
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#' @param neg.log a logical, should a negative logarithm be applied to the posterior copying probabilities to make them distances
#' @param standardize a logical; should the columns be centered and scaled to have unit variance (done after neg.log transform)
#' @param symmetrize a logical; should the final result be symmetried ( for a matrix M, 0.5*(M + t(M)) )
#' @param beta.theta.opts a list; see Details.
#' @return
#'   Matrix of distances.
#'   The \eqn{(j,i)}-th element of of the returned matrix is the inferred
#'   distance \eqn{d_{ji}}{d_(j,i)} between haplotypes \eqn{j} and \eqn{i} at
#'   the current locus.
#'
#' @seealso
#'   \code{\link{PostProbs}} to calculate the posterior marginal probabilities
#'   \eqn{p_{ji}}{p_(j,i)};
#'   \code{\link{Forward}} to propogate a Forward table to a new locus;
#'   \code{\link{Backward}} to propogate a Backward table to a new locus.
#'
#' @examples
#' \dontrun{
#' # To get the distances at, say, locus 10 ...
#' Forward(fwd, pars, 10, nthreads = 8)
#' Backward(bck, pars, 10, nthreads = 8)
#' d <- DistMat(fwd, bck)
#' }
#'
#' @export DistMat
DistMat <- function(fwd, bck, neg.log = TRUE, standardize = FALSE, symmetrize = TRUE, beta.theta.opts = list(bias = 0.5)){

  # Make sure table is not off of the end of the map
  if(fwd$l == 0){"forward table has not been initialized but not propagated to a locus in {1,...,L}."}
  if(bck$l == 2147483647){"backward table has not been initialized but not propagated to a locus in {1,...,L}."}

  if(fwd$pars.sha256 != bck$pars.sha256) {
    stop("parameters used to calculate the forward table and backward table do not match.")
  }

  if(fwd$l > bck$l){stop("fwd$l > bck$l.  The forward table cannot be past the backward table.")}

  if(bck$beta.theta){
    if(fwd$l == bck$l){stop("A forward table cannot be combined with a backward table at the same locus if the Backward Table is in beta.theta space.")}


    if(!is.null(beta.theta.opts$rho1step.fwd) & !is.null(beta.theta.opts$rho1step.bck)){
      # Use custom morgan distances

      if(!is.numeric(beta.theta.opts$rho1step.fwd) || beta.theta.opts$rho1step.fwd<=0 ){stop("rho1step.fwd must be numeric and strictly positive.  To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}
      if(!is.numeric(beta.theta.opts$rho1step.bck) || beta.theta.opts$rho1step.bck<=0 ){stop("rho1step.bck must be numeric and strictly positive.  To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}

      rho.fwd <- beta.theta.opts$rho1step.fwd
      rho.bck <- beta.theta.opts$rho1step.bck

    }else{
      # Use bias

      if(is.null(beta.theta.opts$rho1step.fwd)){stop("rho1step.fwd is missing.")}
      if(is.null(beta.theta.opts$rho1step.bck)){stop("rho1step.bck is missing.")}

      if(!is.numeric(bias) || beta.theta.opts$bias<0 || beta.theta.opts$bias>1 ){stop("bias must be numeric and within [0,1]. To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}

      total.rho <- sum(pars$pars$rho[fwd$l:(bck$l - 1)])

      rho.fwd <- total.rho * bias
      rho.bck <- total.rho * (1-bias)

    }

    tempmat <- ((1 - rho.fwd) * sweep(fwd$alpha,2,fwd$alpha.f,"/") + rho.fwd * pars$pars$Pi) * ((1 - rho.bck) * sweep(bck$beta,2,bck$beta.g,"/") + rho.bck)

  }else{
    if(length(beta.theta.opts)!=0){warning("bck is not in beta.theta space: ignoring beta.theta.opts...")}
    if(bck$l != fwd$l){stop("locus position of the forward table and backward table do not match.")}

    tempmat <- fwd$alpha*bck$beta

  }

  ### Transform

  if(neg.log & standardize){
    tempmat <- -log(tempmat)
    diag(tempmat) <- NA
    tempmat <- scale(tempmat)
    diag(tempmat) <- 0
  }

  if(!neg.log & standardize){
    diag(tempmat) <- NA
    tempmat <- scale(tempmat)
    diag(tempmat) <- 0
  }

  if(neg.log & !standardize){
    tempmat <- sweep(-log(tempmat), MARGIN = 2, STATS = log(colSums(tempmat)), FUN = "+")
    diag(tempmat) <- 0
  }

  if(!neg.log & !standardize){
    tempmat <- sweep(tempmat, MARGIN = 2, STATS = colSums(tempmat), FUN = "/")
  }

  ### Symmetrize

  if(symmetrize){
    tempmat <- 0.5 * (tempmat + t(tempmat))
  }

  class(tempmat) <- c("kalisDistanceMatrix", class(tempmat))

  tempmat
}





#' Plotting function for a kalisDistanceMatrix object
#'
#' Clusters the given distance matrix and generates a heatmap to display it.
#'
#' @param d a kalisDistanceMatrix
#'
#' @return There is nothing returned.
#'
#' @export
plot.kalisDistanceMatrix <- function(x, ...) {
  perm <- fastcluster::hclust(stats::as.dist(x),method="average")$order
  print(lattice::levelplot(x[perm,][,rev(perm)],
                           useRaster = TRUE,
                           col.regions = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,name = "BuPu"))(100),
                           yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxt = "n"))
}
