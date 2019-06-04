#' Calculate recombination parameter rho from more standard genetics parameters
#'
#' Sets up the Li & Stephens HMM transition probabilities rho to be used for a problem.
#'
#' See page 3 in Supplemental Information for the original ChromoPainter paper for motivation behind our parameterization:
#' Lawson, Hellenthal, Myers, and Falush (2012), "Inference of population structure using dense haplotype data", PLoS Genetics, 8 (e1002453).
#'
#' [Insert link to paper for rho maths.]
#'
#' @param morgan.dist a vector specifying the recombination distance between loci in Morgans .
#'   Note element i of this vector should be the distance between loci i and i+1
#'   (not i and i-1), and thus length one less than the haplotype length.  Can be easily obtained by applying \code{diff} to a recombination map "CDF".
#' @param Ne a scalar for the effective population size.
#' @param gamma a scalar power to which the Morgan distances are raised.
#' @param floor if TRUE (default) any transition probabilities below machine precision (1e-16)
#'   will be zeroed out.  If FALSE raw transition probabilities will be preserved.
#'
#' @return A vector of transition probabilities
#'
#' @seealso \code{\link{Parameters}} to use the resulting transition probabilities to construct a kalisParameters object
#'
#' @examples
#' \dontrun{
#' pars <- CalcRho(...)
#' }

#' @export
CalcRho <- function(morgan.dist = 0, Ne = 1, gamma = 1, floor = TRUE) {
  L <- get("L", envir = pkgVars)
  if(anyNA(L)) {
    stop("No haplotypes cached ... cannot determine rho length until cache is loaded with CacheAllHaplotypes().")
  }

  if(!is.vector(floor) || !is.logical(floor) || length(floor) != 1) {
    stop("floor must be TRUE or FALSE.")
  }
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

  # Compute rho
  rho <- c(-(expm1(-Ne*morgan.dist^gamma)), 1)
  if(floor) {
    rho <- ifelse(rho < 1e-16, 0, rho)
  }

  rho
}


#' Create rho from standard parameters
#'
#' Convenience function for calculating kalis recombination (renewal) probabilities from standard Pop. Gen. parameters.
#'
#' Detailed description
#'
#' @param rho recombination probability vector (must be L-1 long), see \code{\link{CalcRho}}
#' @param mu a scalar (for uniform) or vector (for varying) mutation costs.
#' @param Pi leaving the default of uniform copying probabilities is recommended for
#'   computational efficiency.  If desired, a full matrix of background copying
#'   probabilities can be provided, such that the (j,i)-th element is the background
#'   probability that j is copied by i.  Hence, (a) the diagonal must be zero; and (b)
#'   the columns of Pi must sum to 1.
#' @param check.rho if TRUE, a check that rho is within machine precision will be
#'   performed.  If you have created rho using \code{\link{CalcRho}} with \code{floor=TRUE}
#'   then this will be satisfied.
#'
#' @return A \code{kalisParameters} object.
#'
#' @seealso \code{\link{MakeForwardTable}}, \code{\link{MakeForwardTable}} which
#'   construct table objects which internally reference a parameters environment.
#'
#' @examples
#' \dontrun{
#' pars <- Parameters(...)
#' }
#'
#' @export Parameters
Parameters <- function(rho = rep(0, get("L", envir = pkgVars)-1), mu = 1e-8, Pi = NULL, check.rho = TRUE) {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

  if(!is.vector(rho) || length(rho) != L) {
    stop("rho must be a vector of the same length as the sequences in the cache.")
  }
  if(!is.numeric(mu)) {
    stop("rho must be numeric.")
  }
  if(check.rho && any(rho < 1e-16)) {
    stop("some elements of rho are below machine precision.  To disable this check and continue anyway use check.rho=FALSE.")
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
  res$pars$rho <- rho
  res$pars$mu  <- mu
  res$pars$Pi  <- Pi

  # Lock down and checksum
  lockEnvironment(res$pars, bindings = TRUE)
  res$sha256 <- digest(res$pars, algo = "sha256")
  lockEnvironment(res, bindings = TRUE)
  class(res) <- c("kalisParameters", class(res))

  res
}

#' @export
print.kalisParameters <- function(x, ...) {
  if(is.matrix(x$pars$Pi)) {
    Pi <- glue("large matrix, first row: ({glue_collapse(head(x$pars$Pi[1,], 3), ', ')}, ..., {glue_collapse(tail(x$pars$Pi[1,], 3), ', ')})")
  } else {
    Pi <- x$pars$Pi
  }
  if(length(x$pars$mu)>1) {
    mu <- glue("({glue_collapse(head(x$pars$mu, 3), ', ')}, ..., {glue_collapse(tail(x$pars$mu, 3), ', ')})")
  } else {
    mu <- x$pars$mu
  }
  cat(glue("Parameters object with:\n",
           "{'  '}rho   = ({glue_collapse(head(x$pars$rho, 3), ', ')}, ..., {glue_collapse(tail(x$pars$rho, 3), ', ')})\n",
           "{'  '}mu    = {mu}\n",
           "{'  '}Pi    = {Pi}"))
}
