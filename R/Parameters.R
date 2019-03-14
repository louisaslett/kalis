#' Setup genetics parameters
#'
#' Sets up the genetics parameters to be used for a problem.
#'
#' Detailed description
#'
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
#'
#' @return A new parameters environment.
#'
#' @seealso \code{\link{MakeForwardTable}}, \code{\link{MakeForwardTable}} which
#'   construct table objects which internally reference a parameters environment.
#'
#' @examples
#' pars <- Parameters(...)
#'
#' @export Parameters
Parameters <- function(morgan.dist = 0, Ne = 1, gamma = 1, mu = 1e-8, Pi = NULL) {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

  if(!is.numeric(morgan.dist)) {
    stop("morgan.dist must be numeric vector type.")
  }
  if(length(morgan.dist) == 1){
    morgan.dist <- rep(morgan.dist, L-1)
  }
  if(!is.vector(morgan.dist)) {
    stop("morgan.dist must be a vector of recombination distances.")
  }
  if(length(morgan.dist) != L-1) {
    stop("morgan.dist is the wrong length for this problem.")
  }

  if(!is.vector(Ne) || !is.numeric(Ne) || length(Ne) != 1) {
    stop("Ne must be a scalar.")
  }

  if(!is.vector(gamma) || !is.numeric(gamma) || length(gamma) != 1 || gamma < 0) {
    stop("gamma must be a non-negative scalar.")
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

  if(is.null(Pi)) {
    Pi <- 1/(N-1)
  }
  if(is.data.frame(Pi)) {
    stop("Pi must be a matrix or left as NULL to have uniform copying probabilities of 1/(N-1) for a problem with N recipients, not a data frame.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.vector(Pi) && is.numeric(Pi) && length(Pi) == 1 && Pi == 1/(N-1))) {
    stop("Pi must be a matrix or left as NULL to have uniform copying probabilities of 1/(N-1) for a problem with N recipients.")
  }

  # Create the new parameter environment
  res <- new.env(parent = emptyenv())
  res$pars <- new.env(parent = emptyenv())
  res$pars$morgan.dist <- morgan.dist
  res$pars$Ne          <- Ne
  res$pars$gamma       <- gamma
  res$pars$mu          <- mu
  res$pars$Pi          <- Pi

  # Compute rho ... this is the derived parameter we actually care about
  res$pars$rho <- c(1 - exp(-Ne*morgan.dist^gamma), 1)
  res$pars$rho <- ifelse(res$pars$rho < 1e-16, 1e-16, res$pars$rho)

  # Lock down and checksum
  lockEnvironment(res$pars, bindings = TRUE)
  res$sha256 <- digest(res$pars, algo = "sha256")
  lockEnvironment(res, bindings = TRUE)
  class(res) <- c("kalisParameters", class(res))

  res
}

#' @export print.kalisParameters
print.kalisParameters <- function(x, ...) {
  cat(glue("Parameters object with:\n",
           "{'  '}morgan.dist = ({glue_collapse(head(x$pars$morgan.dist, 3), ', ')}, ..., {glue_collapse(tail(x$pars$morgan.dist, 3), ', ')})\n",
           "{'  '}Ne          = {x$pars$Ne}\n",
           "{'  '}gamma       = {x$pars$gamma}\n",
           "{'  '}mu          = {x$pars$mu}\n",
           "{'  '}Pi          = {x$pars$Pi}"))
}
