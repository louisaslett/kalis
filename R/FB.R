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
Forward <- function(fwd, pars, t = fwd$l+1, nthreads = 1) {
  if(!("kalisForwardTable" %in% class(fwd))) {
    stop("The fwd argument is not a valid forward table.")
  }
  if(!("kalisParameters" %in% class(pars))) {
    stop("The pars argument is not a valid parameters object.")
  }
  if(fwd$pars.sha256 != pars$sha256) {
    stop("The forward table provided was created with different parameter values (SHA-256 mismatch).")
  }
  L <- get("hap_size", envir = pkgCache)
  N <- length(get("haps", envir = pkgCache))
  if(!is.vector(t) || !is.numeric(t) || length(t) != 1 || is.na(t)) {
    stop("t must be a scalar.")
  }
  if(t > L) {
    stop(glue("Valid target loci range from 1 to {L} ... cannot move forward to locus {t}."))
  }
  if(fwd$l > t) {
    stop(glue("The forward table provided is for locus position {fwd$l} which is already past requested locus {t}"))
  }
  if(nrow(fwd$alpha) != N || ncol(fwd$alpha) != fwd$to_recipient-fwd$from_recipient+1) {
    stop("Forward table is of the wrong dimensions for this problem.")
  }

  if(is.matrix(pars$pars$Pi)) {
    if(length(pars$pars$mu) == 1) {
      Forward_densePi_scalarmu_cpp(fwd, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, nthreads)
    } else {
      Forward_densePi_densemu_cpp(fwd, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, nthreads)
    }
  } else {
    if(length(pars$pars$mu) == 1) {
      Forward_scalarPi_scalarmu_cpp(fwd, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, nthreads)
    } else {
      Forward_scalarPi_densemu_cpp(fwd, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, nthreads)
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
Backward <- function(bck, pars, t = bck$l-1, nthreads = 1) {
  if(!("kalisBackwardTable" %in% class(bck))) {
    stop("The bck argument is not a valid backward table.")
  }
  if(!("kalisParameters" %in% class(pars))) {
    stop("The pars argument is not a valid parameters object.")
  }
  L <- get("hap_size", envir = pkgCache)
  N <- length(get("haps", envir = pkgCache))
  if(!is.vector(t) || !is.numeric(t) || length(t) != 1 || is.na(t)) {
    stop("t must be a scalar.")
  }
  if(t < 1) {
    stop(glue("Valid target loci range from 1 to {L} ... cannot move backward to locus {t}."))
  }
  if(bck$l < t) {
    stop(glue("The backward table provided is for locus position {bck$l} which is already before requested locus {t}"))
  }
  if(nrow(bck$beta) != N || ncol(bck$beta) != bck$to_recipient-bck$from_recipient+1) {
    stop("Backward table is of the wrong dimensions for this problem.")
  }

  if(is.matrix(pars$pars$Pi)) {
    if(length(pars$pars$mu) == 1) {
      Backward_densePi_scalarmu_cpp(bck, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, nthreads)
    } else {
      Backward_densePi_densemu_cpp(bck, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, nthreads)
    }
  } else {
    if(length(pars$pars$mu) == 1) {
      Backward_scalarPi_scalarmu_cpp(bck, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, nthreads)
    } else {
      Backward_scalarPi_densemu_cpp(bck, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, nthreads)
    }
  }
}
