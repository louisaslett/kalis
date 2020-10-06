#' Propagate a forward table
#'
#' Propagates a \code{kalisForwardTable} to a downstream locus position.
#' The table is updated in-place.
#'
#' \code{Forward} implements the forward algorithm to advance the Li and
#' Stephens rescaled hidden Markov model forward probabilities to a new target
#' locus.
#'
#' @param fwd a \code{kalisForwardTable} object, as returned by
#'   \code{\link{MakeForwardTable}}.
#' @param pars a \code{kalisParameters} object, as returned by
#'   \code{Parameters}.
#' @param t a locus position to move the forward table to.
#'   Must be greater than or equal to current locus position of \code{fwd}.
#'   By default, it simply advances to the next locus.
#' @param nthreads the number of CPU cores to use.
#'   By default no parallelism is used.
#'
#' @return
#' There is nothing returned.
#' For performance reasons, \code{fwd} is updated in-place.
#'
#' @seealso \code{\link{MakeForwardTable}} to generate a forward table;
#'   \code{\link{Backward}} for the analagous backward propagation function;
#'   \code{\link{CopyTable}} to create a copy of table.
#'
#' @examples
#' \dontrun{
#' # Assuming that pars already contains a Parameters object, we can create a
#' # new forward table and propagate to locus 100 using 8-core parallelsim:
#' fwd <- MakeForwardTable()
#' Forward(fwd, pars, 100, nthreads = 8)
#' }
#'
#' @export
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
  L <- get("L", envir = pkgVars)
  N <- get("N", envir = pkgVars)
  if(!is.vector(t) || !is.numeric(t) || length(t) != 1 || is.na(t)) {
    stop("t must be a scalar.")
  }
  if(t < 1 || t > L) {
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
      Forward_densePi_scalarmu_cpp(fwd, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads)
    } else {
      Forward_densePi_densemu_cpp(fwd, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads)
    }
  } else {
    if(length(pars$pars$mu) == 1) {
      Forward_scalarPi_scalarmu_cpp(fwd, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads)
    } else {
      Forward_scalarPi_densemu_cpp(fwd, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads)
    }
  }
}


#' Propagate a backward table
#'
#' Propagates a \code{kalisBackwardTable} to an upstream locus position.
#' The table is updated in-place.
#'
#' \code{Backward} implements the backward algorithm to advance the Li and
#' Stephens rescaled hidden Markov model backward probabilities to a new target
#' locus.
#'
#' @param bck a \code{kalisBackwardTable} object, as returned by
#'   \code{\link{MakeBackwardTable}}.
#' @param pars a \code{kalisParameters} object, as returned by
#'   \code{Parameters}.
#' @param t a locus position to move the backward table to.
#'   Must be less than or equal to current locus position of \code{bck}.
#'   By default, it simply advances to the locus immediately before.
#' @param nthreads the number of CPU cores to use.
#'   By default no parallelism is used.
#' @param beta.theta logical indicating whether the table should be returned in
#'   beta-theta space or in the standard space upon reaching the target locus
#'   \code{t}.
#'
#' @return
#' There is nothing returned.
#' For performance reasons, \code{bck} is updated in-place.
#'
#' @seealso \code{\link{MakeBackwardTable}} to generate a backward table;
#'   \code{\link{Forward}} for the analagous forward propagation function;
#'   \code{\link{CopyTable}} to create a copy of table.
#'
#' @examples
#' \dontrun{
#' # Assuming that pars already contains a Parameters object, we can create a
#' # new backward table and propagate to locus 100 using 8-core parallelsim:
#' bck <- MakeBackwardTable()
#' Backward(bck, pars, 100, nthreads = 8)
#' }
#'
#' @export Backward
Backward <- function(bck, pars, t = bck$l-1, nthreads = 1, beta.theta = FALSE) {
  if(!("kalisBackwardTable" %in% class(bck))) {
    stop("The bck argument is not a valid backward table.")
  }
  if(!("kalisParameters" %in% class(pars))) {
    stop("The pars argument is not a valid parameters object.")
  }
  L <- get("L", envir = pkgVars)
  N <- get("N", envir = pkgVars)
  if(!is.vector(t) || !is.numeric(t) || length(t) != 1 || is.na(t)) {
    stop("t must be a scalar.")
  }
  if(!is.vector(beta.theta) || !is.logical(beta.theta) || length(beta.theta) != 1 || is.na(beta.theta)) {
    stop("beta.theta must be a single logical value.")
  }
  if(t < 1 || t > L) {
    stop(glue("Valid target loci range from 1 to {L} ... cannot move backward to locus {t}."))
  }
  if(bck$l < t) {
    stop(glue("The backward table provided is for locus position {bck$l} which is already before requested locus {t}"))
  }
  if(nrow(bck$beta) != N || ncol(bck$beta) != bck$to_recipient-bck$from_recipient+1) {
    stop("Backward table is of the wrong dimensions for this problem.")
  }
  if(t == bck$l && bck$beta.theta && !beta.theta) {
    stop("Cannot move from beta-theta space to rescaled probability space without moving at least one locus.")
  }

  if(is.matrix(pars$pars$Pi)) {
    if(length(pars$pars$mu) == 1) {
      Backward_densePi_scalarmu_cpp(bck, beta.theta, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads)
    } else {
      Backward_densePi_densemu_cpp(bck, beta.theta, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads)
    }
  } else {
    if(length(pars$pars$mu) == 1) {
      Backward_scalarPi_scalarmu_cpp(bck, beta.theta, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads)
    } else {
      Backward_scalarPi_densemu_cpp(bck, beta.theta, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads)
    }
  }
}
