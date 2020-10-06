#' Calculate recombination parameter
#'
#' A convenient function to calculate the recombination parameter rho from more
#' standard genetics parameters.
#'
#' This convenient function to calculate the recombination parameter rho from
#' other more standard genetics parameters is intended to assist in setting up
#' the Li and Stephens hidden Markov model transition probabilities.
#'
#' See page 3 in Supplemental Information for the original ChromoPainter paper
#' (Lawson et al., 2012) for motivation behind our parameterization, which is as
#' follow:
#'
#' \deqn{\rho = \exp(-N_e \times morgan.dist^\gamma) - 1}{\rho = exp(-Ne * morgan.dist^\gamma) - 1}
#'
#'
#' @param morgan.dist a vector specifying the recombination distance between
#'   loci in Morgans.
#'   Note element i of this vector should be the distance between loci i and i+1
#'   (not i and i-1), and thus length one less than the haplotype length.
#'   This can be easily obtained by applying \code{\link{diff}} to a
#'   recombination map 'CDF'.
#' @param Ne a scalar for the effective population size.
#' @param gamma a scalar power to which the Morgan distances are raised.
#' @param floor if TRUE (default) any transition probabilities below machine
#'   precision (1e-16) will be zeroed out.
#'   If \code{FALSE} raw transition probabilities will be preserved.
#'
#' @return A vector of transition probabilities which can be used at the
#'   \code{rho} argument to the \code{\link{Parameters}} function when creating
#'   a parameter set.
#'
#' @seealso \code{\link{Parameters}} to use the resulting transition
#'   probabilities to construct a \code{kalisParameters} object.
#'
#' @references
#'   Lawson, D. J., Hellenthal, G., Myers, S., & Falush, D. (2012). Inference of
#'   population structure using dense haplotype data. *PLoS genetics*, **8**(1).
#'
#' @examples
#' \dontrun{
#' pars <- CalcRho(...)
#' }
#'
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
  rho <- c(-(expm1(-Ne*morgan.dist^gamma)))
  if(floor) {
    rho <- ifelse(rho < 1e-16, 0, rho)
  }

  rho
}


#' Construct parameter specification
#'
#' Specify a parameter set to be used for a particular Li and Stephens hidden
#' Markov model run.
#'
#' There are 3 parameters which must be specified before the forward or backward
#' equations in the Li and Stephens hidden Markov model can be run.  These are
#' the vector of recombination probabilities, the mutation costs, and the
#' prior copying probabilities.
#'
#' **Recombination probabilities, \code{rho}**
#'
#' This is a vector parameter which must have length \code{L-1}, where \code{L}
#' is the length (number of loci) of the haplotype sequences which have been
#' loaded into the memory cache (using \code{\link{CacheHaplotypes}}).
#' Note that element i of this vector should be the recombination probability
#' between locus \code{i} and \code{i+1}.
#'
#' There is a utility function, \code{\link{CalcRho}}, to assist with creating
#' this recombination vector based on the more standard genetic parameters of
#' Morgan distances, effective population size and a scalar power.
#'
#' By default, the recombination probabilities are set to zero everywhere.
#'
#' **Mutation costs, \code{mu}**
#'
#' The mutation costs may be specified either as uniform across the whole length
#' of the haplotypes (by providing a single scalar value), or may be varying
#' at each locus (by providing a vector of length equal to the length of the
#' haplotype sequences which have been loaded into the memory cache).
#'
#' By default, mutation costs are set to \eqn{10^{-8}}{10^-8}.
#'
#' **Copying probabilities, \code{Pi}**
#'
#' The original Li and Stephens model assumed that each haplotype has an equal
#' prior probability of copying from any other.
#' However, in the spirit of ChromoPainter (Lawson et al., 2012) we allow a
#' matrix of prior copying probabilities.
#'
#' The copying probabilities may be specified as a matrix of size \eqn{N \times N}{N * N}
#' (where \eqn{N} is the number of haplotypes which have been loaded into the
#' memory cache), or for uniform copying probabilities need not be specified
#' (resulting in copying probability \eqn{\frac{1}{N-1}}{1/(N-1)} everywhere).
#' Note that the diagonal must by definition be zero and columns must sum to
#' one.
#'
#' Note that there is a computational cost associated with non-uniform copying
#' probabilities, so absent a strong motivation it is recommended to leave the
#' default of uniform probabilities (NB *do not* specify a uniform matrix which
#' would incur a computation cost too).
#'
#' @param rho recombination probability vector (must be L-1 long).
#'   See \code{\link{CalcRho}} for assistance constructing this from standard
#'   genetics parameters.
#' @param mu a scalar (for uniform) or vector (for varying) mutation costs.
#' @param Pi leaving the default of uniform copying probabilities is recommended
#'   for computational efficiency.
#'   If desired, a full matrix of background copying probabilities can be
#'   provided, such that the (j,i)-th element is the background probability that
#'   j is copied by i.
#'   Hence, (a) the diagonal must be zero; and (b) the columns of Pi must sum to
#'   1.
#' @param check.rho if \code{TRUE}, a check that rho is within machine precision
#'   will be performed.
#'   If you have created rho using \code{\link{CalcRho}} with \code{floor=TRUE}
#'   then this will be satisfied automatically.
#'
#' @return A \code{kalisParameters} object, suitable for use to create the
#'   standard forward and backward recursion tables with
#'   \code{\link{MakeForwardTable}} and \code{\link{MakeBackwardTable}}.
#'   Note you will also need to provide this parameters object when propagating
#'   those tables using either \code{\link{Forward}} or \code{\link{Backward}}.
#'
#' @seealso \code{\link{MakeForwardTable}} and \code{\link{MakeForwardTable}}
#'   which construct table objects which internally reference a parameters
#'   environment.
#'   \code{\link{Forward}} and \code{\link{Backward}} which propagate those
#'   tables according to the Li and Stephens model.
#'
#' @references
#'   Lawson, D. J., Hellenthal, G., Myers, S., & Falush, D. (2012). Inference of
#'   population structure using dense haplotype data. *PLoS genetics*, **8**(1).
#'
#' @examples
#' \dontrun{
#' # To use all the defaults
#' pars <- Parameters()
#'
#' # To use indirect genetics parameters which you have already specified in
#' # morgan.dist, Ne and gamma variables to specify rho, leaving mu and Pi as
#' # defaults
#' pars <- Parameters(CalcRho(morgan.dist, Ne, gamma))
#'
#' # Then use to create a table (eg for forward) ...
#' fwd <- MakeForwardTable(pars)
#' # ... and to propagate it
#' Forward(fwd, pars, 10, nthreads = 8)
#' }
#'
#' @export
Parameters <- function(rho = rep(0, get("L", envir = pkgVars)-1),
                       mu = 1e-8,
                       Pi = NULL,
                       use.speidel = FALSE,
                       check.rho = TRUE) {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

  if(!is.vector(rho) || length(rho) != L-1) {
    stop("rho must be a vector of the length one less than the length of the sequences in the cache.")
  }
  if(!is.numeric(mu)) {
    stop("rho must be numeric.")
  }
  if(check.rho && any(rho < 1e-16 & rho != 0)) {
    stop("some elements of rho are below machine precision but not zero.  To disable this check and continue anyway use check.rho=FALSE.")
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

  if(!test_logical(use.speidel, any.missing = FALSE, len = 1)) {
    stop("use.speidel must be TRUE or FALSE")
  }

  # Create the new parameter environment
  res <- new.env(parent = emptyenv())
  res$pars <- new.env(parent = emptyenv())
  res$pars$rho          <- c(rho, 1)
  res$pars$mu           <- mu
  res$pars$Pi           <- Pi
  res$pars$use.speidel  <- use.speidel

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
