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
#' @param unif.on.underflow a logical; if TRUE, then if all probabilities in a column underflow, then they will be set to \eqn{1/(N-1)} instead of 0
#' @param M an pre-existing matrix into which to write the probabilities, can yield substantial speed up but requires special attention (see Details)
#' @param beta.theta.opts a list; see Details.
#' @param from_recipient offset for distributed problems, see kalis distributed.
#' @param nthreads the number of CPU cores to use. By default no parallelism is used.
#'
#' @return
#'   Matrix of posterior marginal probabilities.
#'   The \eqn{(j,i)}-th element of of the returned matrix is the probability
#'   that \eqn{j} is copied by \eqn{i} at the current locus, \eqn{l}, of the two
#'   tables, given the haplotypes observed (over the whole sequence).
#'
#'
#'   If calculating many posterior probability matrices in succession, providing a pre-existing matrix \code{M} that can be updated in-place can drastically increase speed by eliminating the time needed for memory allocation.
#'   Be warned, since the matrix is updated in-place, if any other variables point to the same memory address, they will also be simultaneously overwritten.  For example, writing
#'   ```
#'   M <- matrix(0,nrow(fwd$alpha),ncol(fwd$alpha))
#'   P <- M
#'   PostProbs(fwd,bck,M=M)
#'   ```
#'   will update M and P simultaneously.
#'
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
PostProbs <- function(fwd, bck, unif.on.underflow = FALSE, M = NULL, beta.theta.opts = NULL, from_recipient = 1, nthreads = 1) {

  rho.list <- input_checks_for_probs_and_dist_mat(fwd,bck,beta.theta.opts)

  # Make M if needed
  if(is.null(M)){ M <- matrix(0,nrow=nrow(fwd$alpha),ncol=ncol(fwd$alpha)) }

  # All clear to calculate distance matrices
  if(bck$beta.theta){
    vector.biproduct <- kalis:::MatAndMulBtwVar(M,fwd,bck,rep(1,nrow(M)),FALSE,TRUE,unif.on.underflow, rho.list$rho.fwd, rho.list$rho.bck, from_recipient,nthreads)
  } else {
    vector.biproduct <- .Call(CCall_MatAndMul, M, fwd, bck, rep(1, nrow(M)), FALSE, TRUE, unif.on.underflow, from_recipient, nthreads)
  }

  invisible(M)
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
#' \deqn{d_{ji} = -log(p_{ji})}{d_(j,i) = - log(p_(j,i)) }
#'
#' where \eqn{p_{ji}}{p_(j,i)} is the posterior marginal probability that
#' \eqn{j} is coped by \eqn{i} at the current locus of the two tables, \eqn{l},
#' given the haplotypes observed (over the whole sequence).
#'
#' By convention, \eqn{d_{ii} = 0}{d_(i,i) = 0} for all \eqn{i}.
#'
#' This function also allows users to calculate distance matrices in between variants and also to calculate matrices that exclude a set of consecutive variants by passing a
#' backward table in beta.theta space.  If in beta.theta space, \code{bck$l} may be greater than but not equal to \code{fwd$l}.  \code{beta.theta.opts} provides is required in this case to
#' set how much of a recombination distance to propagate each matrix before combining them into distances.  See Details below.
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#' @param standardize a logical; should the columns be centered and scaled to have unit variance
#' @param M an pre-existing matrix into which to write the distances, can yield substantial speed up but requires special attention (see Details)
#' @param beta.theta.opts a list; see Details.
#' @param from_recipient offset for distributed problems, see kalis distributed.
#' @param nthreads the number of CPU cores to use. By default no parallelism is used.
#' @return
#'
#'   Matrix of distances.
#'   The \eqn{(j,i)}-th element of of the returned matrix is the inferred
#'   distance \eqn{d_{ji}}{d_(j,i)} to haplotype \eqn{j} from haplotype \eqn{i} at
#'   the current locus.  Each column encodes the output of an independent HMM: in column \eqn{i}, haplotype \eqn{i} is taken as the observed recipient haplotype and painted as a mosaic of the other \eqn{N-1} haplotypes.
#'   Hence, the distances are asymmetric.
#'
#'   If calculating many distance matrices in succession, providing a pre-existing matrix \code{M} that can be updated in-place can drastically increase speed by eliminating the time needed for memory allocation.
#'   Be warned, since the matrix is updated in-place, if any other variables point to the same memory address, they will also be simultaneously overwritten.  For example, writing
#'   ```
#'   M <- matrix(0,nrow(fwd$alpha),ncol(fwd$alpha))
#'   P <- M
#'   DistMat(fwd,bck,M=M)
#'   ```
#'   will update M and P simultaneously.
#'
#'   When provided, \code{M} must have dimensions matching that of \code{fwd$alpha}.  Typically, that is simply \eqn{N} x \eqn{N} for \eqn{N} haplotypes.
#'   However, if kalis is being run distributed, \code{M} will be a \eqn{N} x \eqn{R} matrix where \eqn{R} is the number of recipient haplotypes on the given machine.
#'
#'
#'   In order to obtain distance matrices between variants \code{fwd$l} and \code{bck$l}, then \code{bck} must be in \code{beta.theta} space.  This allows the forward and backward tables to be transitioning both tables to some genomic position between \code{fwd$l} and \code{bck$l}.
#'   The precise recombination distance by which each table is propagated can be controlled by passing optional arguments in a list via \code{beta.theta.opts}.
#'   The recombination distances used can be specified in one of two ways.
#'
#'   1. Manually.  In this case, \code{beta.theta.opts} is a list containing two named elements:
#'      -  \code{"rho.fwd"} \eqn{\in (0,1)} specifies the transition probability \eqn{rho} for propagating the forward table
#'      -  \code{"rho.bck"} \eqn{\in (0,1)} specifies the transition probability \eqn{rho} for propagating the backward table.
#'
#'   2. Implicitly.  In this case, \code{beta.theta.opts} is a list containing two named elements:
#'      - \code{"pars"}: a \code{kalisParameters} object that implicitly defines the recombination distance \eqn{\rho^\star} between \code{fwd$l} and \code{bck$l}
#'      - \code{"bias"} \eqn{\in (0,1)}.  The forward table is propagated a distance of \code{bias}\eqn{\rho^\star} and the backward table is propagated a distance of \code{(1-bias)}\eqn{\rho^\star}.
#'
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
DistMat <- function(fwd, bck, standardize = FALSE, M = NULL, beta.theta.opts = NULL, from_recipient = 1, nthreads = 1){

  rho.list <- input_checks_for_probs_and_dist_mat(fwd,bck,beta.theta.opts)

  # Make M if needed
  if(is.null(M)){ M <- matrix(0,nrow=nrow(fwd$alpha),ncol=ncol(fwd$alpha)) }


  # All clear to calculate distance matrices
  if(bck$beta.theta){
    vector.biproduct <- kalis:::MatAndMulBtwVar(M,fwd,bck,rep(1,nrow(M)),standardize, FALSE, FALSE, rho.list$rho.fwd, rho.list$rho.bck, from_recipient,nthreads)
  } else {
    vector.biproduct <- .Call(CCall_MatAndMul, M, fwd, bck, rep(1, nrow(M)), standardize, FALSE, FALSE, from_recipient, nthreads)
  }

  invisible(M)
}



input_checks_for_probs_and_dist_mat <-  function(fwd,bck,beta.theta.opts){

  # RUN GENERAL CHECKS
  if(fwd$l == 0){stop("forward table has not been initialized but not propagated to a locus in {1,...,L}.")}
  if(bck$l == 2147483647){stop("backward table has not been initialized but not propagated to a locus in {1,...,L}.")}

  if(fwd$pars.sha256 != bck$pars.sha256) {
    stop("parameters used to calculate the forward table and backward table do not match.")
  }

  if(fwd$l > bck$l){stop("fwd$l > bck$l.  The forward table cannot be past the backward table.")}


  # RUN BTW LOCI and AT LOCI specific checks
  if(bck$beta.theta){

    if(fwd$l == bck$l){stop("A forward table cannot be combined with a backward table at the same locus if the backward table is in beta.theta space.")}

    if(!((!is.null(beta.theta.opts$rho.fwd) & !is.null(beta.theta.opts$rho.bck)) | (!is.null(beta.theta.opts$pars) & !is.null(beta.theta.opts$bias)))){
      stop("beta.theta.opts must be a named list containing either pars and bias OR rho.fwd and rho.bck.")
    }

    if(!is.null(beta.theta.opts$rho.fwd) & !is.null(beta.theta.opts$rho.bck)){
      # Use custom morgan distances

      if(!is.numeric(beta.theta.opts$rho.fwd) || beta.theta.opts$rho.fwd<=0 ){stop("rho.fwd must be numeric and strictly positive.  To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}
      if(!is.numeric(beta.theta.opts$rho.bck) || beta.theta.opts$rho.bck<=0 ){stop("rho.bck must be numeric and strictly positive.  To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}

      rho.fwd <- beta.theta.opts$rho.fwd
      rho.bck <- beta.theta.opts$rho.bck

    }else{
      # Use bias and pars

      if(!inherits(beta.theta.opts$pars,"kalisParameters")){stop("beta.theta.opts$pars must be kalisParameters object.")}

      if(!is.numeric(beta.theta.opts$bias) || beta.theta.opts$bias<=0 || beta.theta.opts$bias>=1 ){stop("bias must be numeric and within (0,1). To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}

      total.rho <- sum(beta.theta.opts$pars$pars$rho[fwd$l:(bck$l - 1)])

      rho.fwd <- total.rho * beta.theta.opts$bias
      rho.bck <- total.rho * (1-beta.theta.opts$bias)

    }

    return(list("rho.fwd" = rho.fwd, "rho.bck" = rho.bck))

  }else{

    if(bck$l != fwd$l){stop("locus position of the forward table and backward table do not match.")}
    return(NULL)
  }
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
