#' Calculate recombination probabilities
#'
#' Calculate the recombination probabilities, rho, from a given recombination map for the haplotype data currently in the kalis cache.
#'
#' This is a utility function to calculate the recombination probabilities rho (the Li and Stephens hidden Markov model transition probabilities) from a recombination map/distances.
#'
#' **NOTE:** the corresponding haplotype data *must* have already been inserted into the kalis cache by a call to [CacheHaplotypes()], since this function performs checks to confirm the dimensionality matches.
#'
#' See page 3 in Supplemental Information for the original ChromoPainter paper (Lawson et al., 2012) for motivation behind our parameterisation, which is as follows:
#'
#' \deqn{\rho = 1 - \exp(-s \times cM^\gamma)}{\rho = 1 - exp(-s * cM^\gamma)}
#'
#' For a complete description, see the main kalis paper, Aslett and Christ (2024).
#'
#' @references
#' Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation of the Li & Stephens model for local ancestry inference in R", *BMC Bioinformatics*, **25**(1). Available at: \doi{10.1186/s12859-024-05688-8}.
#'
#' Lawson, D.J., Hellenthal, G., Myers, S. and Falush, D. (2012) "Inference of population structure using dense haplotype data", *PLoS genetics*, **8**(1). Available at: \doi{10.1371/journal.pgen.1002453}.
#'
#' @param cM a vector specifying the recombination distance between variants in centimorgans.
#'   Note element i of this vector should be the distance between variants `i` and `i+1` (not `i` and `i-1`), and thus length one less than the number of variants.
#'   This can be easily obtained by applying [diff()] to a recombination map 'CDF'.
#'   By default, recombination probabilities are zero, meaning that all variants are perfectly linked.
#' @param s a scalar multiplier on the recombination map (related to effective population size).
#' @param gamma a scalar power to which the Morgan distances are raised.
#' @param floor if `TRUE` (default) any recombination probabilities below machine precision (`1e-16`) will be zeroed out.
#'   If `FALSE` raw recombination probabilities will be preserved.
#'
#' @return A vector of recombination probabilities which can be used as the `rho` argument to the [Parameters()] function when creating a parameter set.
#'
#' @seealso [Parameters()] to use the resulting recombination probabilities to construct a `kalisParameters` object.
#'
#' @examples
#' # Load the mini example data and recombination map from the package built-in #' # dataset
#' data("SmallHaps")
#' data("SmallMap")
#'
#' # Load haplotypes into cache
#' CacheHaplotypes(SmallHaps)
#'
#' # Compute the recombination probabilities, leaving s and gamma as default
#' rho <- CalcRho(diff(SmallMap))
#' rho
#'
#' @export
CalcRho <- function(cM = 0, s = 1, gamma = 1, floor = TRUE) {
  L <- get("L", envir = pkgVars)
  if(anyNA(L)) {
    stop("No haplotypes cached ... cannot determine rho length until cache is loaded with CacheAllHaplotypes().")
  }

  if(!is.vector(floor) || !is.logical(floor) || length(floor) != 1) {
    stop("floor must be TRUE or FALSE.")
  }
  if(!is.numeric(cM)) {
    stop("cM must be numeric vector type.")
  }
  if(length(cM) == 1){
    cM <- rep(cM, L-1)
  }
  if(!is.vector(cM)) {
    stop("cM must be a vector of recombination distances.")
  }
  if(length(cM) != L-1) {
    stop("cM is the wrong length for this problem.")
  }
  morgans <- cM/100.0

  if(!test_number(s, lower = 0, finite = TRUE)) {
    stop("s must be a positive scalar.")
  }

  if(!is.vector(gamma) || !is.numeric(gamma) || length(gamma) != 1 || gamma < 0) {
    stop("gamma must be a non-negative scalar.")
  }

  # Compute rho
  rho <- c(-(expm1(-s*morgans^gamma)))
  if(floor) {
    rho <- ifelse(rho < 1e-16, 0, rho)
  }

  rho
}


#' Define a set of Li and Stephens parameters
#'
#' Specify a parameter set to be used for a particular Li and Stephens hidden Markov model run.
#'
#' There are 3 parameters which must be specified before the forward or backward equations in the Li and Stephens hidden Markov model can be run.
#' These are the vector of recombination probabilities, the mutation probabilities, and the prior copying probabilities.
#'
#' **NOTE:** the corresponding haplotype data *must* have already been inserted into the kalis cache by a call to [CacheHaplotypes()], since this function performs checks to confirm the dimensionality matches.
#'
#' **Recombination probabilities, `rho`**
#'
#' This is a vector parameter which must have length \eqn{L-1}, where \eqn{L} is the number of variants that have been loaded into the kalis memory cache (using [CacheHaplotypes()]).
#' Note that element `i` of this vector should be the recombination probability between variants `i` and `i+1`.
#'
#' There is a utility function, [CalcRho()], to assist with creating these recombination probabilities from a recombination map.
#'
#' By default, the recombination probabilities are set to zero everywhere.
#'
#' **Mutation probabilities, `mu`**
#'
#' The mutation probabilities may be specified either as uniform across all variants (by providing a single scalar value), or may be varying at each variant (by providing a vector of length equal to the number of variants, \eqn{L}, which have been loaded into the kalis memory cache).
#'
#' By default, mutation probabilities are set to \eqn{10^{-8}}{10^-8}.
#'
#' **Copying probabilities, `Pi`**
#'
#' The original Li and Stephens model assumed that each haplotype has an equal prior probability of copying from any other.
#' However, in the spirit of ChromoPainter (Lawson et al., 2012) we allow a matrix of prior copying probabilities.
#'
#' The copying probabilities may be specified as a standard R matrix of size \eqn{N \times N}{N x N} (where \eqn{N} is the number of haplotypes which have been loaded into the kalis memory cache).
#' The element at row `j`, column `i` corresponds to the prior (background) probability that haplotype `i` copies from haplotype `j`.
#' Note that the diagonal must by definition be zero and columns must sum to one.
#' Alternatively, for uniform copying probabilities, this argument need not be specified (resulting in copying probability \eqn{\frac{1}{N-1}}{1/(N-1)} everywhere).
#'
#' Note that there is a computational cost associated with non-uniform copying probabilities, so it is recommended to leave the default of uniform probabilities when appropriate (**Note:** *do not* specify a uniform matrix when uniform probabilities are intended, since this would end up incurring the computational cost of non-uniform probabilities).
#'
#' @references
#' Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation of the Li & Stephens model for local ancestry inference in R", *BMC Bioinformatics*, **25**(1). Available at: \doi{10.1186/s12859-024-05688-8}.
#'
#' Lawson, D.J., Hellenthal, G., Myers, S.R. and Falush, D. (2012) "Inference of population structure using dense haplotype data", *PLoS Genetics*, **8**(1). Available at: \doi{10.1371/journal.pgen.1002453}.
#'
#' Speidel, L., Forest, M., Shi, S. and Myers, S.R. (2019) "A method for genome-wide genealogy estimation for thousands of samples", *Nature Genetics*, **51**, p. 1321-1329. Available at: \doi{10.1038/s41588-019-0484-x}.
#'
#' @param rho recombination probability vector (must be \eqn{L-1} long).
#'   See [CalcRho()] for assistance constructing this from a recombination
#'   map/distances.
#' @param mu a scalar (for uniform) or vector (for varying) mutation probabilities.
#' @param Pi leaving the default of uniform copying probabilities is recommended
#'   for computational efficiency.
#'   If desired, a full matrix of background copying probabilities can be provided, such that the `[j,i]`-th element is the background probability that `i` copies from `j`.
#'   Hence, (a) the diagonal must be zero; and (b) the columns of `Pi` must sum to 1.
#'   Note: each column corresponds to an independent Li and Stephens hidden Markov model.
#' @param use.speidel a logical, if `TRUE`, use the asymmetric mutation model used by RELATE (Speidel et al., 2019).
#'   **WARNING:** this model assumes that the cached haplotypes have an ancestral/derived encoding -- zeros denote ancestral variant carriers and ones denote derived variant carriers. Defaults to `FALSE`.
#' @param check.rho if `TRUE`, a check that rho is within machine precision will be performed.
#'   If you have created rho using [CalcRho()] with `floor = TRUE` then this will be satisfied automatically.
#'   This can be important to ensure that the convex combination `rho` and `(1-rho)` normalises to 1 within machine precision, but an advanced user may wish to override this requirement.
#'
#' @return A `kalisParameters` object, suitable for use to create the standard forward and backward recursion tables with [MakeForwardTable()] and [MakeBackwardTable()].
#'   Note you will also need to provide this parameters object when propagating those tables using either [Forward()] or [Backward()].
#'
#' @seealso [MakeForwardTable()] and [MakeBackwardTable()] which construct table objects which internally reference a parameters environment;
#'   [Forward()] and [Backward()] which propagate those tables according to the Li and Stephens model.
#'
#' @examples
#' # Load the mini example data and recombination map from the package built-in #' # dataset
#' data("SmallHaps")
#' data("SmallMap")
#'
#' # Load haplotypes into cache
#' CacheHaplotypes(SmallHaps)
#'
#' # To use all the defaults
#' pars <- Parameters()
#' pars
#'
#' # Or perhaps use SmallMap to compute the recombination probabilities, leaving
#' # s and gamma as default
#' rho <- CalcRho(diff(SmallMap))
#' pars <- Parameters(rho)
#' pars
#'
#' # Then use these parameters to create a table (eg for forward) ...
#' fwd <- MakeForwardTable(pars)
#' # ... and to propagate it
#' Forward(fwd, pars, 10)
#' fwd
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

#' @importFrom utils tail
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
           "{'  '}Pi    = {Pi}"), "\n")
}
