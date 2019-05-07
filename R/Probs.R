#' Calculate Posterior Marginal Probabilities from a Forward Table Object and a Backward Table Object
#'
#' Provides an easy function for calculating the posterior marginal probabilities.
#'
#' The forward and backward tables must be at the same locus in order for them to be combined to yield the posterior marginal probabilities at locus l.
#' The (j,i)-th element of of the returned matrix is the probability that j is copied by i at the current locus of the two tables, l, given the haplotypes observed (from locus 1 to L).
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
#' \code{\link{DistMat}} to generate calculate \eqn{d_{ji}} distances directly
#' \code{\link{MakeForwardTable}} to generate Forward table;
#' \code{\link{Forward}} to propogate a Backward table to a new locus;
#' \code{\link{MakeBackwardTable}} to generate a Backward table;
#' \code{\link{Backward}} to propogate a Backward table to a new locus.
#'
#' @examples
#' \dontrun{
#' PostProbs(fwd,bck)
#' }
#'
#' @export PostProbs
PostProbs <- function(fwd, bck, log = FALSE) {
  if(fwd$l != bck$l) {
    warning("Computing dist matrix but locus position of the forward table and backward table do not match.")
  }
  if(fwd$pars.sha256 != bck$pars.sha256) {
    warning("Computing dist matrix but parameters used to calculate the forward table and backward table do not match.")
  }
  tempmat <- fwd$alpha*bck$beta
  tempmat <- sweep(log(tempmat), MARGIN = 2, STATS = log(colSums(tempmat)), FUN = "-")
  if(log) {
    return(tempmat)
  } else {
    return(exp(tempmat))
  }
}



#' Calculate a Distance Matrix from a Forward Table and a Backward Table
#'
#' Calculates a local distance matrix from a \code{kalisForwardTable} and a \code{kalisBackwardTable} at the same locus.
#'
#' The forward and backward tables must be at the same locus l.
#' The (j,i)-th element of of the returned matrix is the inferred distance d_(j,i) between haplotypes j and i at the current locus of the two tables, l, given the haplotypes observed (from locus 1 to L).
#'
#' d_(j,i) = -( log(p_(j,i)) + log(p_(i,j)) ) / 2 where p_(j,i) is the posterior marginal probability that j is coped by i at the current locus of the two tables, l, given the haplotypes observed (from locus 1 to L).
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
#' \code{\link{PostProbs}} to calculate the posterior marginal probabilities p_(j,i);
#' \code{\link{MakeForwardTable}} to generate Forward table;
#' \code{\link{Forward}} to propogate a Backward table to a new locus;
#' \code{\link{MakeBackwardTable}} to generate a Backward table;
#' \code{\link{Backward}} to propogate a Backward table to a new locus.
#'
#' @examples
#' \dontrun{
#' DistMat(fwd, bak)
#' }
#'
#' @export DistMat
DistMat <- function(fwd, bck) {
  if(fwd$l != bck$l) {
    warning("Computing dist matrix but locus position of the forward table and backward table do not match.")
  }
  if(fwd$pars.sha256 != bck$pars.sha256) {
    warning("Computing dist matrix but parameters used to calculate the forward table and backward table do not match.")
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
