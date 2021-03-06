#' Propagate a forward table
#'
#' Propagates a \code{kalisForwardTable} to a downstream variant.
#' The table is updated in-place.
#'
#' \code{Forward} implements the forward algorithm to advance the Li and
#' Stephens rescaled hidden Markov model forward probabilities to a new target
#' variant.  Naturally, this can only propagate a table to variants downstream
#' of its current position.
#'
#' For mathematical details please see Section 2 of the kalis paper (TODO: ref).
#' Note that the precise formulation of the forward equation is determined by
#' whether the flag \code{use.spiedel} is set in the parameters provided in
#' \code{pars}.
#'
#' @param fwd a \code{kalisForwardTable} object, as returned by
#'   \code{\link{MakeForwardTable}}.
#' @param pars a \code{kalisParameters} object, as returned by
#'   \code{\link{Parameters}}.
#' @param t a variant to move the forward table to.
#'   Must be greater than or equal to current variant of \code{fwd}.
#'   By default, it simply advances to the next variant downstream.
#' @param nthreads the number of CPU cores to use.
#'   By default uses the \code{parallel} package to detect the number of physical cores.
#'
#' @return
#' There is nothing returned.
#' For performance reasons, \code{fwd} is updated in-place.
#'
#' @seealso \code{\link{MakeForwardTable}} to generate a forward table;
#'   \code{\link{Backward}} for the analogous backward propagation function;
#'   \code{\link{CopyTable}} to create a copy of table.
#'
#' @examples
#' \dontrun{
#' # Load the toy haplotype example and set toy parameters
#' CacheHaplotypes(SmallHaplotypes)
#' pars <- Parameters(rho = CalcRho(cM = SmallMap))
#'
#' # Create a forward table for the hidden Markov model incorporating all
#' # recipient and donor haplotypes
#' fwd <- MakeForwardTable(pars)
#'
#' # Create a forward table for the hidden Markov model incorporating only
#' # recipient haplotypes 100 to 200 (inclusive) and all donor haplotypes.
#' fwd <- MakeForwardTable(pars, 100, 200)
#'
#' # This table is uninitialised, but ready to pass to the Forward function
#' # which will trigger initialisation and propagation.
#' # For example, initialise and propagate forward to the 10th variant:
#' Forward(fwd, pars, 10, nthreads = 8)
#' print(fwd)
#' }
#'
#' @export
Forward <- function(fwd,
                    pars,
                    t = fwd$l+1,
                    nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient-fwd$from_recipient+1)) {
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
  if(!test_integerish(t, lower = 1, upper = L, any.missing = FALSE, len = 1)) {
    stop(glue("t must be a scalar integer in the range from 1 to {L}."))
  }
  t <- as.integer(t)
  if(fwd$l <= L && fwd$l > t) {
    stop(glue("The forward table provided is for locus position {fwd$l} which is already past requested locus {t}"))
  }
  if(nrow(fwd$alpha) != N || ncol(fwd$alpha) != fwd$to_recipient-fwd$from_recipient+1) {
    stop("Forward table is of the wrong dimensions for this problem.")
  }

  if(identical(nthreads, "R")) {
    warning("Warning: using gold master R implementation.")
    return(invisible(Forward.GM(fwd, pars, t)))
  }
  maxthreads <- parallel::detectCores()
  if(!test_integerish(nthreads, lower = 1, upper = maxthreads, any.missing = FALSE, len = 1) &&
     !test_integerish(nthreads, lower = 0, upper = maxthreads-1, any.missing = FALSE, min.len = 2, unique = TRUE)) {
    stop(glue("The nthreads argument must either be a single scalar (between 1 and your number of logical cores (which is {maxthreads})) indicating the number of threads, or a vector indicating the core number to run each thread on (so in total length(nthreads) threads) are run."))
  }
  if((length(nthreads) > 1 && length(nthreads) > fwd$to_recipient-fwd$from_recipient+1) || (length(nthreads) == 1 && nthreads > fwd$to_recipient-fwd$from_recipient+1)) {
    stop(glue("Cannot launch more threads than there are recipients in the forward table (here to_recipient to from_recipient covers {fwd$to_recipient-fwd$from_recipient+1} recipients)."))
  }
  nthreads <- as.integer(nthreads)

  invisible(.Call(CCall_Forward, fwd, FALSE, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads))
}


#' Propagate a backward table
#'
#' Propagates a \code{kalisBackwardTable} to an upstream variant.
#' The table is updated in-place.
#'
#' \code{Backward} implements the backward algorithm to propagate the Li and
#' Stephens rescaled hidden Markov model backward probabilities to a new target
#' variant.  Naturally, this can only propagate a table to variants upstream
#' of its current position.
#'
#' For mathematical details please see Section 2 of the kalis paper (TODO: ref).
#' Note that the precise formulation of the backward equation is determined by
#' whether the flag \code{use.spiedel} is set in the parameters provided in
#' \code{pars}.
#'
#' **Beta-theta space**
#'
#' The rescaled HMM backward probabilities incorporate
#' all of the haplotype relatedness information downstream from but NOT including the
#' target variant, as is standard in the definition of HMM backward probabilities -- we refer to this as "beta space."
#' In contrast, the rescaled forward probabilities, naturally incorporate all of the
#' haplotype relatedness information upstream from AND including the target variant.
#' Setting \code{beta.theta = TRUE} incorporates the target variant (analogous to the rescaled forward probabilities) --
#' we refer to this as "beta-theta space".
#'
#' A backward table in beta-theta space (with \code{beta.theta = TRUE}) can be propagated
#' to an upstream variant without incorporating that variant, thereby moving to beta space (\code{beta.theta = FALSE}), and vice versa.
#' However, while a backward table in beta space (\code{beta.theta = FALSE})
#' can be updated to incorporate the current variant, a backward table
#' that is already in beta-theta space can not move to beta space
#' without changing variants -- that would involve "forgetting" the current variant (see Examples).
#'
#'
#' @param bck a \code{kalisBackwardTable} object, as returned by
#'   \code{\link{MakeBackwardTable}}.
#' @param pars a \code{kalisParameters} object, as returned by
#'   \code{\link{Parameters}}.
#' @param t a variant to move the backward table to.
#'   Must be less than or equal to current variant of \code{bck}.
#'   By default, it simply advances to the variant immediately upstream.
#' @param nthreads the number of CPU cores to use.
#'   By default uses the \code{parallel} package to detect the number of physical cores.
#' @param beta.theta logical indicating whether the table should be returned in
#'   beta-theta space or in the standard space upon reaching the target variant
#'   \code{t}.  See the details section.
#'
#' @return
#' There is nothing returned.
#' For performance reasons, \code{bck} is updated in-place.
#'
#' @seealso \code{\link{MakeBackwardTable}} to generate a backward table;
#'   \code{\link{Forward}} for the analogous forward propagation function;
#'   \code{\link{CopyTable}} to create a copy of table.
#'
#' @examples
#' \dontrun{
#' # Load the toy haplotype example and set toy parameters
#' CacheHaplotypes(SmallHaplotypes)
#' pars <- Parameters(rho = CalcRho(cM = SmallMap))
#'
#' # Create a backward table for the hidden Markov model incorporating all
#' # recipient and donor haplotypes
#' bck <- MakeBackwardTable(pars)
#'
#' # Create a backward table for the hidden Markov model incorporating only
#' # recipient haplotypes 100 to 200 (inclusive) and all donor haplotypes.
#' bck <- MakeBackwardTable(pars, 100, 200)
#'
#' # This table is uninitialised, but ready to pass to the Backward function
#' # which will trigger initialisation and propagation from the last variant.
#' # For example, initialise and propagate backward to the 10th variant:
#' Backward(bck, pars, 10, nthreads = 8)
#' print(bck)
#'
#' #### Beta-theta space example ####
#' # Now moving to beta-theta space
#' Backward(bck, pars, 8, beta.theta=TRUE)
#'
#' # Now moving to beta space (can be done with or without propagating to a new variant)
#' Backward(bck, pars, 8, beta.theta=FALSE)
#'
#' # Attempting to move from beta space to beta-theta space without propagating
#' # is not possible (see Details).
#' Backward(bck, pars, 8, beta.theta=TRUE)
#' }
#'
#' @export Backward
Backward <- function(bck,
                     pars,
                     t = bck$l-1,
                     nthreads = min(parallel::detectCores(logical = FALSE), bck$to_recipient-bck$from_recipient+1),
                     beta.theta = FALSE) {
  if(!("kalisBackwardTable" %in% class(bck))) {
    stop("The bck argument is not a valid backward table.")
  }
  if(!("kalisParameters" %in% class(pars))) {
    stop("The pars argument is not a valid parameters object.")
  }
  if(bck$pars.sha256 != pars$sha256) {
    stop("The backward table provided was created with different parameter values (SHA-256 mismatch).")
  }
  L <- get("L", envir = pkgVars)
  N <- get("N", envir = pkgVars)
  if(!test_integerish(t, lower = 1, upper = L, any.missing = FALSE, len = 1)) {
    stop(glue("t must be a scalar integer in the range from 1 to {L}."))
  }
  t <- as.integer(t)
  if(bck$l < t) {
    stop(glue("The backward table provided is for locus position {bck$l} which is already before requested locus {t}"))
  }
  if(!is.vector(beta.theta) || !is.logical(beta.theta) || length(beta.theta) != 1 || is.na(beta.theta)) {
    stop("beta.theta must be a single logical value.")
  }
  if(nrow(bck$beta) != N || ncol(bck$beta) != bck$to_recipient-bck$from_recipient+1) {
    stop("Backward table is of the wrong dimensions for this problem.")
  }
  if(t == bck$l && bck$beta.theta && !beta.theta) {
    stop("Cannot move from beta-theta space to beta space without moving at least one locus.")
  }

  if(identical(nthreads, "R")) {
    warning("Warning: using gold master R implementation.")
    return(invisible(Backward.GM(bck, pars, t, beta.theta)))
  }
  maxthreads <- parallel::detectCores()
  if(!test_integerish(nthreads, lower = 1, upper = maxthreads, any.missing = FALSE, len = 1) &&
     !test_integerish(nthreads, lower = 0, upper = maxthreads-1, any.missing = FALSE, min.len = 2, unique = TRUE)) {
    stop(glue("The nthreads argument must either be a single scalar (between 1 and your number of logical cores (which is {maxthreads})) indicating the number of threads, or a vector indicating the core number to run each thread on (so in total length(nthreads) threads) are run."))
  }
  if((length(nthreads) > 1 && length(nthreads) > bck$to_recipient-bck$from_recipient+1) || (length(nthreads) == 1 && nthreads > bck$to_recipient-bck$from_recipient+1)) {
    stop(glue("Cannot launch more threads than there are recipients in the forward table (here to_recipient to from_recipient covers {bck$to_recipient-bck$from_recipient+1} recipients)."))
  }
  nthreads <- as.integer(nthreads)

  invisible(.Call(CCall_Backward, bck, beta.theta, t, pars$pars$Pi, pars$pars$mu, pars$pars$rho, pars$pars$use.speidel, nthreads))
}
