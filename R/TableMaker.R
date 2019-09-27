#' Create a \code{kalisForwardTable}
#'
#' Allocates the memory for and initializes a forward table.
#'
#' \code{MakeForwardTable} returns a \code{kalisForwardTable} object appropriate for
#' a given set of haplotypes (that must have already been cached by \code{CacheAllHaplotypes}) and
#' a given set of HMM parameters specified by \code{pars}.  The returned \code{kalisForwardTable} is
#' initialized at locus 0 and is ready to be propagated to a given target locus with the function \code{\link{Forward}}.
#'
#' Since there is an independent hidden Markov model run for each recipient
#' haplotype, it is possible to create a partial forward table object which
#' corresponds to a subset of recipients using the \code{from_recipient} and
#' \code{to_recipient} arguments.
#'
#' @param pars a \code{kalisParameters} object specifying the genetics
#'   parameters to be associated with this forward table.  These parameters can
#'   be set up by using the \code{\link{Parameters}} function.
#' @param from_recipient first recipient haplotype if creating a partial forward
#'   table.  By default includes from the first recipient haplotype.
#' @param to_recipient last recipient haplotype if creating a partial forward
#'   table.  By default includes to the last recipient haplotype.
#'
#' @return A specialized list of class \code{kalisForwardTable}.  For a given \code{kalisForwardTable}, \code{fwd},
#' \code{fwd$l} denotes the current locus position of \code{fwd}.  \code{fwd$alpha} is a matrix of rescaled forward probabilities under the Li \& Stephens HMM.
#' Each column of \code{fwd$alpha} corresponds to an independent HMM such that \eqn{\alpha^l_{ji}}
#' is proportional to the probability that haplotype \eqn{j} is copied by haplotype \eqn{i} at locus \eqn{l} and
#' observing haplotype \eqn{i} from locus 1 up through locus \eqn{l}.  \code{fwd$alpha.f}
#' is a vector containing scaling constants needed to continue propagating the HMM (please see kalis paper for details).
#'
#' \code{kalisForwardTable} also carries with it a checksum key for the parameters \code{pars} it was provided.
#' If a user attempts to interact a \code{kalisForwardTable} with a \code{kalisBackwardTable} with mismatched
#' parameters, an error will be thrown.
#'
#'
#' @seealso \code{\link{Forward}} to propagate the newly created \code{kalisForwardTable}.
#'
#' @examples
#' # Examples
#' \dontrun{
#' # Create a forward table for the hidden Markov model incorporating all
#' # recipient and donor haplotypes
#' fwd <- MakeForwardTable()
#'
#' # Create a forward table for the hidden Markov model incorporating only
#' # recipient haplotypes 100 to 200 (inclusive) and all donor haplotypes.
#' fwd <- MakeForwardTable(100, 200)
#'
#' # This table is uninitialised, but ready to pass to the Forward function
#' # which will trigger initialisation and propagation from the first locus.
#' # For example, initialise and propagate forward to locus 10:
#' Forward(fwd, 10, morgan.dist, Ne, gamma, mu, nthreads = 8)
#' }
#'
#' @aliases print.kalisForwardTable
#' @export MakeForwardTable
MakeForwardTable <- function(pars, from_recipient = 1, to_recipient = Inf) {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

  if(!("kalisParameters" %in% class(pars))) {
    if("rstudioapi" %in% utils::installed.packages()[, "Package"]) {
      rstudioapi::sendToConsole("?Parameters")
    }
    stop("The pars argument is not a valid parameters object.  See Parameters() function for how to create it.")
  }

  if(from_recipient>to_recipient) {
    stop("from_recipient must be smaller than to_recipient.")
  }
  if(from_recipient < 1) {
    from_recipient <- 1L
  }
  if(to_recipient > N) {
    to_recipient <- N
  }
  if(as.integer(from_recipient) != from_recipient) {
    stop("from_recipient must be an integer.")
  }
  if(as.integer(to_recipient) != to_recipient) {
    stop("to_recipient must be an integer.")
  }
  from_recipient <- as.integer(from_recipient)
  to_recipient <- as.integer(to_recipient)
  delN <- to_recipient-from_recipient+1

  fwd <- list()

  # Define core table, duplicating where relevant to ensure unique to this forward table
  fwd$alpha          <- duplicate(matrix(0, N, delN))
  fwd$alpha.f        <- duplicate(rep(0, delN))
  fwd$l              <- duplicate(c(0L))
  fwd$from_recipient <- duplicate(from_recipient)
  fwd$to_recipient   <- duplicate(to_recipient)
  fwd$pars.sha256    <- duplicate(pars$sha256)

  class(fwd) <- c("kalisForwardTable", class(fwd))

  fwd
}

#' @export
print.kalisForwardTable <- function(x, ...) {
  if(!("kalisForwardTable" %in% class(x)))
    stop("Not a kalisForwardTable object")

  d <- dim(x$alpha)

  if(d[1]==d[2]) {
    cat(glue("Full Forward Table object for {d[1]} haplotypes."), "\n")
  } else {
    cat(glue("Partial Forward Table object for {d[1]} haplotypes."), "\n")
    cat(glue("  Recipients {x$from_recipient} to {x$to_recipient}"), "\n")
  }

  if(x$l == 0) {
    cat("  Newly created table, currently uninitialised to any locus (ready for Forward function next).\n")
  } else {
    cat(glue("  Current locus = {x$l}"), "\n")
  }
  cat("  Memory consumed: ", ceiling(utils::object.size(x)/1e6)/1e3, "GB.\n")
}

#' Create a \code{kalisBackwardTable}
#'
#' Allocates the memory for and initializes a backward table.
#'
#' \code{MakeBackwardTable} returns a \code{kalisBackwardTable} object appropriate for
#' a given set of haplotypes (that must have already been cached by \code{CacheAllHaplotypes}) and
#' a given set of HMM parameters specified by \code{pars}.  The returned \code{kalisBackwardTable} is
#' initialized at the end of the cached haplotypes (technically \code{bck$l}=2,147,483,647 for computational reasons)
#' and is ready to be propagated to a given target locus with the function \code{\link{Backward}}.
#'
#' Since there is an independent hidden Markov model run for each recipient
#' haplotype, it is possible to create a partial forward table object which
#' corresponds to a subset of recipients using the \code{from_recipient} and
#' \code{to_recipient} arguments.
#'
#' @param pars a \code{kalisParameters} environment specifying the genetics
#'   parameters to be associated with this forward table.  These parameters can
#'   be set up by using the \code{\link{Parameters}} function.
#' @param from_recipient first recipient haplotype if creating a partial forward
#'   table.  By default includes from the first recipient haplotype.
#' @param to_recipient last recipient haplotype if creating a partial forward
#'   table.  By default includes to the last recipient haplotype.
#'
#' @return A specialized list of class \code{kalisBackwardTable}.  For a given \code{kalisBackwardTable}, \code{bck},
#'  \code{bck$l} denotes the current locus position of \code{fwd}.  \code{bck$beta} is a matrix of rescaled backward probabilities under the Li \& Stephens HMM.
#' Each column of \code{bck$beta} corresponds to an independent HMM such that \eqn{\beta^\ell_{ji}}
#' is proportional to the probability of observing haplotype \eqn{i} from locus \eqn{l+1} up through locus \eqn{L} given that haplotype \eqn{j}
#' is copied by haplotype \eqn{i} at locus \eqn{\ell}.  \code{bck$beta.g} and \code{bck$beta.g2}
#' are both vectors containing scaling constants needed to continue propagating the HMM (please see kalis paper for details).
#'
#' \code{kalisBackwardTable} also carries with it a checksum key for the parameters \code{pars} it was provided.
#' If a user attempts to interact a \code{kalisBackwardTable} with a \code{kalisForwardTable} with mismatched
#' parameters, an error will be thrown.
#'
#
#' @seealso \code{\link{Backward}} to propagate the newly created \code{kalisBackwardTable}.
#'
#' @examples
#' \dontrun{
#' MakeBackwardTable(...)
#' }
#'
#' @aliases print.kalisBackwardTable
#' @export MakeBackwardTable
MakeBackwardTable <- function(pars, from_recipient = 1, to_recipient = Inf) {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

  if(!("kalisParameters" %in% class(pars))) {
    if("rstudioapi" %in% utils::installed.packages()[, "Package"]) {
      rstudioapi::sendToConsole("?Parameters")
    }
    stop("The pars argument is not a valid parameters object.  See Parameters() function for how to create it.")
  }

  if(from_recipient>to_recipient) {
    stop("from_recipient must be smaller than to_recipient.")
  }
  if(from_recipient < 1) {
    from_recipient <- 1L
  }
  if(to_recipient > N) {
    to_recipient <- N
  }
  if(as.integer(from_recipient) != from_recipient) {
    stop("from_recipient must be an integer.")
  }
  if(as.integer(to_recipient) != to_recipient) {
    stop("to_recipient must be an integer.")
  }
  from_recipient <- as.integer(from_recipient)
  to_recipient <- as.integer(to_recipient)
  delN <- to_recipient-from_recipient+1

  bck <- list()

  # Define core table, duplicating where relevant to ensure unique to this forward table
  bck$beta           <- duplicate(matrix(0, N, delN))
  bck$beta.g         <- duplicate(rep(0, delN))
  bck$beta.g2        <- duplicate(rep(0, delN))
  bck$l              <- duplicate(c(2147483647L))
  bck$from_recipient <- duplicate(from_recipient)
  bck$to_recipient   <- duplicate(to_recipient)
  bck$pars.sha256    <- duplicate(pars$sha256)

  class(bck) <- c("kalisBackwardTable", class(bck))

  bck
}

#' @export
print.kalisBackwardTable <- function(x, ...) {
  if(!("kalisBackwardTable" %in% class(x)))
    stop("Not a kalisBackwardTable object")

  d <- dim(x$beta)

  if(d[1]==d[2]) {
    cat(glue("Full Backward Table object for {d[1]} haplotypes."), "\n")
  } else {
    cat(glue("Partial Backward Table object for {d[1]} haplotypes."), "\n")
    cat(glue("  Recipients {x$from_recipient} to {x$to_recipient}"), "\n")
  }

  if(x$l == 2147483647) {
    cat("  Newly created table, currently uninitialised to any locus (ready for Backward function next).\n")
  } else {
    cat(glue("  Current locus = {x$l}"), "\n")
  }
  cat("  Memory consumed: ", ceiling(utils::object.size(x)/1e6)/1e3, "GB.\n")
}
