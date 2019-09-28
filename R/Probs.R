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
#' Calculates a local distance matrix.
#'
#' This computes a local distance matrix based on the forward and backward
#' tables at a certain locus.
#' The forward and backward tables provided must be at the same locus \eqn{l}.
#' The \eqn{(j,i)}-th element of of the returned matrix is the inferred distance
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
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#'
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
DistMat <- function(fwd, bck) {
  if(fwd$l != bck$l) {
    stop("locus position of the forward table and backward table do not match.")
  }
  if(fwd$pars.sha256 != bck$pars.sha256) {
    stop("parameters used to calculate the forward table and backward table do not match.")
  }

  tempmat <- fwd$alpha*bck$beta
  tempmat <- sweep(-log(tempmat), MARGIN = 2, STATS = log(colSums(tempmat)), FUN = "+")
  diag(tempmat) <- 0

  d <- (tempmat + t(tempmat))/2

  class(d) <- c("kalisDistanceMatrix", class(d))

  d
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
