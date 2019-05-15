#' Propagate a Forward Table
#'
#' Propagates a \code{kalisForwardTable} to a downstream locus position (nothing is returned: the table is updated in-place).
#'
#' Forward implements the forward algorithm to advance the rescaled HMM forward probabilities stored in
#' a \code{kalisForwardTable} \code{fwd} to a new target locus \code{t}.  Note that \code{t} must be greater than the current locus, \code{fwd$l}.
#'
#' @param fwd \code{kalisForwardTable} as returned by \code{\link{MakeForwardTable}}
#' @param pars a \code{kalisParameters} object returned by \code{Parameters}
#' @param t a locus position to move the forward table to.  Must be greater than
#'   or equal to \code{fwd$l}.
#' @param nthreads the number of CPU cores to use
#'
#' @return There is nothing returned.  For performance reasons, \code{fwd} is updated in-place.
#'
#' @seealso \code{\link{MakeForwardTable}} to generate a forward table;
#'   \code{\link{Backward}} for the analagous backward propagation function.
#'
#' @examples
#' \dontrun{
#' fwd <- MakeForwardTable()
#' Forward(fwd, 100, Pi, mu, rho)
#' }
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
  L <- get("L", envir = pkgVars)
  N <- get("N", envir = pkgVars)
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


#' Propagate a Backward Table
#'
#' Propagates a \code{kalisBackwardTable} to an upstream locus position (nothing is returned: the table is updated in-place).
#'
#' Backward implements the backward algorithm to advance the rescaled HMM backward probabilities stored in
#' \code{bck} to a new target locus \code{t}.  Note that \code{t} must be less than or equal to the current locus, \code{bck$l}.
#'
#' @param bck \code{kalisBackwardTable} as returned by \code{\link{MakeBackwardTable}}
#' @param pars a \code{kalisParameters} object returned by \code{Parameters}
#' @param t a target locus position to move the backward table to.  Must be less than or equal to the current locus, \code{bck$l}.
#' @param nthreads the number of CPU cores to use
#'
#' @return There is nothing returned.  For performance reasons, \code{bck} is updated in-place.
#'
#' @seealso \code{\link{MakeBackwardTable}} to generate backward table;
#'   \code{\link{Forward}} for the analagous forward propagation function.
#'
#' @examples
#' \dontrun{
#' bck <- MakeBackwardTable()
#' Backward(bck, pars, 2500)
#' }
#'
#' @export Backward
Backward <- function(bck, pars, t = bck$l-1, nthreads = 1) {
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
