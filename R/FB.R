#' Propagate an HMM forward table
#'
#' Takes a forward table and propagates it in-place to a later locus position.
#'
#' Detailed description
#'
#' @param fwd a forward table as returned by \code{\link{MakeForwardTable}}
#' @param t a locus position to move the forward table to.  Must be greater than
#'   or equal to locus position of table provided in \code{fwd}.
#' @param Pi a matrix of background copying probabilities.  Can also provide a
#'   scalar value for uniform background copying probability.
#' @param mu a vector of ...
#' @param rho a vector of ...
#' @param nthreads the number of CPU cores on which to run
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
Forward <- function(fwd, t, Pi, mu, rho, nthreads = 1) {
  L <- get("seq_size", envir = pkgCache)
  N <- length(get("seqs", envir = pkgCache))
  if(fwd$l > t) {
    stop("The forward table provided is for locus position ", l, " which is already past requested locus ", t)
  }
  if(nrow(fwd$alpha) != N || ncol(fwd$alpha) != fwd$to_recipient-fwd$from_recipient+1) {
    stop("Forward table is of the wrong dimensions for this problem.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.atomic(Pi) && length(Pi) == 1L && !is.character(Pi) && Im(Pi)==0)) {
    stop("If Pi is not a full matrix of background copying probabilities then it must be a single scalar for uniform background copying probability.")
  }
  if(!is.vector(mu) || length(mu) != L) {
    stop("mu is of the wrong type/length.")
  }
  if(!is.vector(rho) || length(rho) != L) {
    stop("rho is of the wrong type/length.")
  }

  if(is.matrix(Pi)) {
    Forward_densePi_cpp(fwd, t, Pi, mu, rho, nthreads)
  } else {
    Forward_scalarPi_cpp(fwd, t, Pi, mu, rho, nthreads)
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
#' @param Pi a matrix of background copying probabilities.  Can also provide a
#'   scalar value for uniform background copying probability.
#' @param mu a vector of ...
#' @param rho a vector of ...
#' @param nthreads the number of CPU cores on which to run
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
Backward <- function(bck, t, Pi, mu, rho, nthreads = 1) {
  L <- get("seq_size", envir = pkgCache)
  N <- length(get("seqs", envir = pkgCache))
  if(bck$l < t) {
    stop("The backward table provided is for locus position ", l, " which is already before requested locus ", t)
  }
  if(nrow(bck$beta) != N || ncol(bck$beta) != bck$to_recipient-bck$from_recipient+1) {
    stop("Forward table is of the wrong dimensions for this problem.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.atomic(Pi) && length(Pi) == 1L && !is.character(Pi) && Im(Pi)==0)) {
    stop("If Pi is not a full matrix of background copying probabilities then it must be a single scalar for uniform background copying probability.")
  }
  if(!is.vector(mu) || length(mu) != L) {
    stop("mu is of the wrong type/length.")
  }
  if(!is.vector(rho) || length(rho) != L) {
    stop("rho is of the wrong type/length.")
  }

  if(is.matrix(Pi)) {
    Backward_densePi_cpp(bck, t, Pi, mu, rho, nthreads)
  } else {
    Backward_scalarPi_cpp(bck, t, Pi, mu, rho, nthreads)
  }
}





