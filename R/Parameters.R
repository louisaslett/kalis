#' Calculate recombination parameter rho from more standard genetics parameters
#'
#' Sets up the Li & Stephens HMM transition probabilities rho to be used for a problem.
#'
#' Detailed description
#'
#' @param morgan.dist a vector specifying the recombination distance between loci in Morgans .
#'   Note element i of this vector should be the distance between loci i and i+1
#'   (not i and i-1), and thus length one less than the haplotype length.  Can be easily obtained by applying \code{diff} to a recombination map "CDF".
#' @param Ne a scalar for the effective population size.
#' @param gamma a scalar power to which the Morgan distances are raised.
#' @param floor.rho a logical, if TRUE (default), then all rho that are initially calculated to be less than 1e-16 are set to zero.  If FALSE, all rho initially calculated to be less than 1e-16 are set to 1e-16.
#'
#' @return A vector of transition probabilities
#'
#' @seealso \code{\link{Parameters}} to use the resulting transition probabilities to construct a kalisParameters object
#'
#' @examples
#' \dontrun{
<<<<<<< HEAD
#' pars <- CalcRho(...)
#' }
#'
#' @export CalcRho
CalcRho <- function(morgan.dist = 0, Ne = 1, gamma = 1, floor.rho = TRUE) {
=======
#' pars <- Parameters(CalcRho(Morgan.dist, 1, 1), rep(1e-08, L), Pi)
#' }
#'
#' @export
Parameters <- function(rho = 1e-16, mu = 1e-8, Pi = NULL) {
>>>>>>> e87a9ec840220b404718e1da1159bc2150b5ea5f
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

  if(!is.numeric(rho)) {
    stop("rho must be numeric vector type.")
  }
  if(!is.vector(rho)) {
    stop("rho must be a vector of recombination distances.")
  }
  if(length(rho) == 1){
    rho <- rep(rho, L)
  }
  if(length(rho) != L) {
    stop("rho is the wrong length for this problem.")
  }


  # Compute rho
  rho <- c(1 - exp(-Ne*morgan.dist^gamma), 1)

  if(floor.rho){
    rho <- ifelse(rho < 1e-16, 0, rho)
  }else{
    rho <- ifelse(rho < 1e-16, 1e-16, rho)
  }

  rho
}



#' Setup kalis HMM parameters
#'
#' Sets up the genetics parameters to be used for a problem.
#'
#' Detailed description
#'
#' @param rho recombination probability vector (must be L-1 long), see \code{CalcRho}
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
#' \dontrun{
#' pars <- Parameters(...)
#' }
#'
#' @export Parameters
Parameters <- function(rho = rep(0,get("L", envir = pkgVars)-1), mu = 1e-8, Pi = NULL) {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

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
CalcRho <- function(morgan.dist = 0, Ne = 1, gamma = 1, threshold = 1e-16) {
  L <- get("L", envir = pkgVars)
  if(anyNA(L)) {
    stop("No haplotypes cached ... cannot determine rho length until cache is loaded with CacheAllHaplotypes().")
  }

  if(threshold < 0) {
    stop("threshold must be non-negative.")
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

  # Compute rho ... this is the derived parameter we actually care about
  rho <- c(1 - exp(-Ne*morgan.dist^gamma), 1)
  rho <- ifelse(rho < threshold, threshold, rho)
  rho
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
