#' Auto-tune the parameters for HMM
#'
#' Provides auto-tuning on a given chromosome for the \code{Ne}, \code{mu} and
#' \code{gamma} parameters of the HMM.
#'
#' Detailed description
#'
#' Also reference the vignette
#'
#' @param t the number of loci to auto-tune against (which will be evently spaced
#'   according to the inverse cumulative recombination map); or a vector of specific
#'   loci to use for auto-tuning.
#' @param cache a forward table cache as returned by \code{\link{CreateForwardTableCache}}.
#'   Note this cache will be overwritten.
#' @param morgan.dist a vector of recombination distances between loci, in Morgans.
#'   Note element i of this vector should be the distance between loci i and i+1
#'   (not i and i-1), and thus length one less than the sequence length.
#' @param Pi leaving the default of uniform copying probabilities is recommended for
#'   computational efficiency.  If desired, a full matrix of background copying
#'   probabilities can be provided, such that the (i,j)-th element is the background
#'   probability that i copies j.  Hence, (a) the diagonal must be zero; and (b)
#'   the columns of Pi must sum to 1.
#' @param nthreads the number of CPU cores on which to run.
#'
#' @return The optimal parameters for the HMM as a list.
#'
#' @seealso \code{\link{Forward}}, \code{\link{ForwardUsingTableCache}}, and
#'   \code{\link{Backward}} for functions which consume the parameters which
#'   can be auto-tuned by this function.
#'
#' @examples
#'
AutoTune <- function(t, cache, morgan.dist, Pi = 1/(nrow(fwd$alpha)-1), nthreads = 1) {
  if(!is.vector(t)) {
    stop("t must be either a vector or a scalar.")
  }
  if(!is.numeric(t)) {
    stop("t must be numeric.")
  }
  if(length(t) < 1 || length(t) > L) {
    stop("t is the wrong length for this problem.")
  }


  if(length(t) == 1) {
    t <-
  }

}

InvRecombMap <- (morgan.dist) {

}
