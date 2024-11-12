#' Propagate a forward table
#'
#' Propagates a `kalisForwardTable` to a downstream variant.
#' The table is updated in-place.
#'
#' `Forward` implements the forward algorithm to advance the Li and Stephens rescaled hidden Markov model forward probabilities to a new target variant.
#' Naturally, this can only propagate a table to variants downstream of its current position.
#'
#' For mathematical details please see Section 2 of the kalis paper (https://doi.org/10.1186/s12859-024-05688-8).
#' Note that the precise formulation of the forward equation is determined by whether the flag `use.spiedel` is set in the parameters provided in `pars`.
#'
#' @param fwd a `kalisForwardTable` object, as returned by
#'   [MakeForwardTable()].
#' @param pars a `kalisParameters` object, as returned by
#'   [Parameters()].
#' @param t a variant to move the forward table to.
#'   Must be greater than or equal to current variant of `fwd`.
#'   By default, it advances to the next variant downstream (or if an uninitialised table, to the first variant).
#' @param nthreads if a scalar, the number of CPU cores to use.
#'   If a vector, launch as many threads as the length of the vector and attempt to pin the threads to those CPU cores (requires system to support thread affinity).
#'   By default uses the `parallel` package to detect the number of physical cores.
#'
#' @return
#' There is nothing returned.
#'
#' **NOTE:** for performance reasons, `fwd` is updated in-place.
#'
#' @seealso [MakeForwardTable()] to generate a forward table;
#'   [Backward()] for the analogous backward propagation function;
#'   [CopyTable()] to create a copy of table.
#'
#' @examples
#' # Load the toy haplotype example and set toy parameters
#' data("SmallHaps")
#' data("SmallMap")
#'
#' CacheHaplotypes(SmallHaps)
#'
#' rho <- CalcRho(diff(SmallMap))
#' pars <- Parameters(rho)
#'
#' # Create the forward table we want to propagate
#' fwd <- MakeForwardTable(pars)
#'
#' # Calling Forward on this uninitialised table moves it to the first variant
#' Forward(fwd, pars)
#' fwd
#'
#' # And again moves it to the next variant (etcetera)
#' Forward(fwd, pars)
#' fwd
#'
#' # Or zoom to a particular variant
#' Forward(fwd, pars, 50)
#' fwd
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
  if(t == 2147483647L + 1) {
    t <- 1L
  }
  if(!test_integerish(t, lower = 1, upper = L, any.missing = FALSE, len = 1)) {
    stop(glue("t must be a scalar integer in the range from 1 to {L}."))
  }
  t <- as.integer(t)
  if(fwd$l <= L && fwd$l > t) {
    stop(glue("The forward table provided is for variant position {fwd$l} which is already past requested variant {t}"))
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
#' Propagates a `kalisBackwardTable` to an upstream variant.
#' The table is updated in-place.
#'
#' `Backward` implements the backward algorithm to propagate the Li and Stephens rescaled hidden Markov model backward probabilities to a new target
#' variant.
#' Naturally, this can only propagate a table to variants upstream of its current position.
#'
#' For mathematical details please see Section 2 of the kalis paper (https://doi.org/10.1186/s12859-024-05688-8).
#' Note that the precise formulation of the backward equation is determined by whether the flag `use.spiedel` is set in the parameters provided in `pars`.
#'
#' **Beta-theta space**
#'
#' The rescaled HMM backward probabilities incorporate all of the haplotype relatedness information downstream from but NOT including the target variant, as is standard in the definition of HMM backward probabilities -- we refer to this as "beta space", or "rescaled probability space."
#' In contrast, the rescaled forward probabilities, naturally incorporate all of the haplotype relatedness information upstream from AND including the target variant.
#' Setting `beta.theta = TRUE` incorporates the target variant (analogous to the rescaled forward probabilities) -- we refer to this as "beta-theta space."
#'
#' A backward table in beta-theta space (with `beta.theta = TRUE`) can be propagated to an upstream variant without incorporating that variant, thereby moving to beta space (`beta.theta = FALSE`), and vice versa.
#' However, while a backward table in beta space (`beta.theta = FALSE`) can be updated to incorporate the current variant, a backward table that is already in beta-theta space can not move to beta space without changing variants -- that would involve "forgetting" the current variant (see Examples).
#'
#'
#' @param bck a `kalisBackwardTable` object, as returned by
#'   [MakeBackwardTable()].
#' @param pars a `kalisParameters` object, as returned by
#'   [Parameters()].
#' @param t a variant to move the backward table to.
#'   Must be less than or equal to current variant of `bck`.
#'   By default, it advances to the next variant upstream (or if an uninitialised table, to the last variant).
#' @param nthreads if a scalar, the number of CPU cores to use.
#'   If a vector, launch as many threads as the length of the vector and attempt to pin the threads to those CPU cores (requires system to support thread affinity).
#'   By default uses the `parallel` package to detect the number of physical cores.
#' @param beta.theta logical indicating whether the table should be returned in beta-theta space or in the standard space upon reaching the target variant `t`.
#'   See the Details section.
#'
#' @return
#' There is nothing returned.
#'
#' **NOTE:** for performance reasons, `bck` is updated in-place.
#'
#' @seealso [MakeBackwardTable()] to generate a backward table;
#'   [Forward()] for the analogous forward propagation function;
#'   [CopyTable()] to create a copy of table.
#'
#' @examples
#'
#' # Load the toy haplotype example and set toy parameters
#' data("SmallHaps")
#' data("SmallMap")
#'
#' CacheHaplotypes(SmallHaps)
#'
#' rho <- CalcRho(diff(SmallMap))
#' pars <- Parameters(rho)
#'
#' # Create the backward table we want to propagate
#' bck <- MakeBackwardTable(pars)
#'
#' # Calling Backward on this uninitialised table moves it to the last variant
#' Backward(bck, pars)
#' bck
#'
#' # And again moves it to the next variant (etcetera)
#' Backward(bck, pars)
#' bck
#'
#' # Or zoom to a particular variant
#' Backward(bck, pars, 150)
#' bck
#'
#' # Now moving to variant 125 AND specifying "beta space" (though this is the
#' # default, just being very clear)
#' Backward(bck, pars, 125, beta.theta = FALSE)
#' bck
#'
#' # Now just moving to "beta-theta space" (can be done with or without
#' # propagating to a new variant)
#' Backward(bck, pars, 125, beta.theta = TRUE)
#' bck
#'
#' # However, attempting to move from "beta-theta space" back to "beta space"
#' # without propagating is not possible (see Details).
#' # The following will give an error (hence wrapped in try())
#' try(Backward(bck, pars, 125, beta.theta = FALSE))
#' bck
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
  if(t == 2147483647L - 1) {
    t <- L
  }
  if(!test_integerish(t, lower = 1, upper = L, any.missing = FALSE, len = 1)) {
    stop(glue("t must be a scalar integer in the range from 1 to {L}."))
  }
  t <- as.integer(t)
  if(bck$l < t) {
    stop(glue("The backward table provided is for variant position {bck$l} which is already before requested variant {t}"))
  }
  if(!is.vector(beta.theta) || !is.logical(beta.theta) || length(beta.theta) != 1 || is.na(beta.theta)) {
    stop("beta.theta must be a single logical value.")
  }
  if(nrow(bck$beta) != N || ncol(bck$beta) != bck$to_recipient-bck$from_recipient+1) {
    stop("Backward table is of the wrong dimensions for this problem.")
  }
  if(t == bck$l && bck$beta.theta && !beta.theta) {
    stop("Cannot move from beta-theta space to beta space without moving at least one variant.")
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
