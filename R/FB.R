#' Propagate an HMM Forward Table
#'
#' Takes a forward table and propagates it in-place to a later locus position.
#'
#' Forward implements the forward algorithm to advance the lag-scaled forward probabilities stored in
#' a forward table fwd to a new target locus t.  Note that t must be greater than the current fwd locus.
#' The standard forward probabilities
#' mention that each column is an independent HMM and highlight that
#' forward table (and Pi) are therefore to be viewed column-wise.
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#' @param t a locus position to move the forward table to.  Must be greater than
#'   or equal to locus position of table provided in \code{fwd}.
#' @param morgan.dist a vector of recombination distances between loci, in Morgans.
#'   Note element i of this vector should be the distance between loci i and i+1
#'   (not i and i-1), and thus length one less than the haplotype length.
#' @param Ne a scalar for the effective population size.  Can be autotuned, see ...
#' @param gamma a scalar power to which the Morgan distances are raised.  Can be
#'   autotuned, see ...
#' @param mu a scalar (for uniform) or vector (for varying) mutation costs.
#' @param Pi leaving the default of uniform copying probabilities is recommended for
#'   computational efficiency.  If desired, a full matrix of background copying
#'   probabilities can be provided, such that the (i,j)-th element is the background
#'   probability that j copies i.  Hence, (a) the diagonal must be zero; and (b)
#'   the columns of Pi must sum to 1.
#' @param nthreads the number of CPU cores on which to run.
#'
#' @return There is nothing returned.  For performance reasons, the forward
#'   table which was passed in is updated in-place.
#'
#' @seealso \code{\link{MakeForwardTable}} to generate forward table;
#'   \code{\link{Backward}} for analagous backward induction function.
#'
#' @examples
#' fwd <- MakeForwardTable()
#' Forward(fwd, 100, Pi, mu, rho)
#'
#' @export Forward
Forward <- function(fwd, t, morgan.dist, Ne, gamma, mu, Pi = 1/(nrow(fwd$alpha)-1), nthreads = 1) {
  L <- get("hap_size", envir = pkgCache)
  N <- length(get("haps", envir = pkgCache))
  if(fwd$l > t) {
    stop("The forward table provided is for locus position ", fwd$l, " which is already past requested locus ", t)
  }
  if(nrow(fwd$alpha) != N || ncol(fwd$alpha) != fwd$to_recipient-fwd$from_recipient+1) {
    stop("Forward table is of the wrong dimensions for this problem.")
  }
  if(!is.vector(morgan.dist)) {
    stop("morgan.dist must be a vector of recombination distances.")
  }
  if(!is.numeric(morgan.dist)) {
    stop("morgan.dist must be numeric vector type.")
  }
  if(length(morgan.dist) != L-1) {
    stop("morgan.dist is the wrong length for this problem.")
  }
  if(!is.vector(Ne) || !is.numeric(Ne) || length(Ne) != 1) {
    stop("Ne must be a scalar.")
  }
  if(!is.vector(gamma) || !is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
    stop("gamma must be a positive scalar.")
  }
  if(!is.vector(mu)) {
    stop("mu must be either a vector or a scalar.")
  }
  if(!is.numeric(mu)) {
    stop("mu must be numeric.")
  }
  if(length(mu) != 1 && length(mu) != L) {
    stop("mu is the wrong length for this problem.")
  }
  if(is.data.frame(Pi)) {
    stop("Pi must be a matrix or scalar, not a data frame.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.vector(Pi) && is.numeric(Pi) && length(Pi) == 1 && Pi == 1/(nrow(fwd$alpha)-1))) {
    stop("Pi can only be set to a matrix, or omitted to have uniform copying probabilities of 1/(N-1) for a problem with N recipients.")
  }

  rho <- c(1-exp(-Ne*morgan.dist^gamma), 1)
  rho <- ifelse(rho<1e-16, 1e-16, rho)

  if(is.matrix(Pi)) {
    if(length(mu) == 1) {
      Forward_densePi_scalarmu_cpp(fwd, t, Pi, mu, rho, nthreads)
    } else {
      Forward_densePi_densemu_cpp(fwd, t, Pi, mu, rho, nthreads)
    }
  } else {
    if(length(mu) == 1) {
      Forward_scalarPi_scalarmu_cpp(fwd, t, Pi, mu, rho, nthreads)
    } else {
      Forward_scalarPi_densemu_cpp(fwd, t, Pi, mu, rho, nthreads)
    }
  }
}

#' Propagate an HMM backward table
#'
#' Takes a backward table and propagates it in-place to an earlier locus
#' position.
#'
#' Detailed description
#'
#' @param bck a backward table as returned by \code{\link{MakeBackwardTable}}
#' @param t a locus position to move the backward table to.  Must be less than
#'   or equal to locus position of table provided in \code{bck}.
#' @param morgan.dist a vector of recombination distances between loci, in Morgans.
#'   Note element i of this vector should be the distance between loci i and i+1
#'   (not i and i-1), and thus length one less than the haplotype length.
#' @param Ne a scalar for the effective population size.  Can be autotuned, see ...
#' @param gamma a scalar power to which the Morgan distances are raised.  Can be
#'   autotuned, see ...
#' @param mu a scalar (for uniform) or vector (for varying) mutation costs.
#' @param Pi leaving the default of uniform copying probabilities is recommended for
#'   computational efficiency.  If desired, a full matrix of background copying
#'   probabilities can be provided, such that the (i,j)-th element is the background
#'   probability that j copies i.  Hence, (a) the diagonal must be zero; and (b)
#'   the columns of Pi must sum to 1.
#' @param nthreads the number of CPU cores on which to run.
#'
#' @return There is nothing returned.  For performance reasons, the backward
#'   table which was passed in is updated in-place.
#'
#' @seealso \code{\link{MakeBackwardTable}} to generate backward table;
#'   \code{\link{Forward}} for analagous forward induction function.
#'
#' @examples
#' bck <- MakeBackwardTable()
#' Backward(bck, 100, Pi, mu, rho)
#'
#' @export Backward
Backward <- function(bck, t, morgan.dist, Ne, gamma, mu, Pi = 1/(nrow(bck$beta)-1), nthreads = 1) {
  L <- get("hap_size", envir = pkgCache)
  N <- length(get("haps", envir = pkgCache))
  if(bck$l < t) {
    stop("The backward table provided is for locus position ", bck$l, " which is already before requested locus ", t)
  }
  if(nrow(bck$beta) != N || ncol(bck$beta) != bck$to_recipient-bck$from_recipient+1) {
    stop("Forward table is of the wrong dimensions for this problem.")
  }
  if(!is.vector(morgan.dist)) {
    stop("morgan.dist must be a vector of recombination distances.")
  }
  if(!is.numeric(morgan.dist)) {
    stop("morgan.dist must be numeric vector type.")
  }
  if(length(morgan.dist) != L-1) {
    stop("morgan.dist is the wrong length for this problem.")
  }
  if(!is.vector(Ne) || !is.numeric(Ne) || length(Ne) != 1) {
    stop("Ne must be a scalar.")
  }
  if(!is.vector(gamma) || !is.numeric(gamma) || length(gamma) != 1) {
    stop("gamma must be a scalar.")
  }
  if(!is.vector(mu)) {
    stop("mu must be either a vector or a scalar.")
  }
  if(!is.numeric(mu)) {
    stop("mu must be numeric.")
  }
  if(length(mu) != 1 && length(mu) != L) {
    stop("mu is the wrong length for this problem.")
  }
  if(is.data.frame(Pi)) {
    stop("Pi must be a matrix or scalar, not a data frame.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.vector(Pi) && is.numeric(Pi) && length(Pi) == 1 && Pi == 1/(nrow(bck$beta)-1))) {
    stop("Pi can only be set to a matrix, or omitted to have uniform copying probabilities of 1/(N-1) for a problem with N recipients.")
  }

  rho <- c(1-exp(-Ne*morgan.dist^gamma), 1)
  rho <- ifelse(rho<1e-16, 1e-16, rho)

  if(is.matrix(Pi)) {
    if(length(mu) == 1) {
      Backward_densePi_scalarmu_cpp(bck, t, Pi, mu, rho, nthreads)
    } else {
      Backward_densePi_densemu_cpp(bck, t, Pi, mu, rho, nthreads)
    }
  } else {
    if(length(mu) == 1) {
      Backward_scalarPi_scalarmu_cpp(bck, t, Pi, mu, rho, nthreads)
    } else {
      Backward_scalarPi_densemu_cpp(bck, t, Pi, mu, rho, nthreads)
    }
  }
}
