#' Create a `kalisForwardTable`
#'
#' Allocates the memory for and initializes a forward table.
#'
#' `MakeForwardTable` initializes a `kalisForwardTable` object that will store the rescaled forward probabilities for each recipient haplotype.
#' Note: since all other haplotypes loaded by [CacheHaplotypes()] into the kalis cache are taken as potential donor haplotypes, this table will be of size (total haplotypes in cache)x(number of recipients).
#' Also note that this table object is linked to a given set of Li and Stephens hidden Markov model parameters, created by [Parameters()], to ensure consistency for any given HMM run.
#'
#' The returned `kalisForwardTable` is ready to be propagated to a given target variant with the function [Forward()].
#'
#' Since each column corresponds to an independent Li and Stephens hidden Markov model (ie for each recipient), it is possible to create a partial forward table object which corresponds to a subset of recipients using the `from_recipient` and `to_recipient` arguments.
#'
#' @param pars a `kalisParameters` object specifying the genetics parameters to be associated with this forward table.
#'   These parameters can be set up by using the [Parameters()] function.
#' @param from_recipient first recipient haplotype included if creating a partial forward table.
#'   By default all are included from the first recipient haplotype.
#'   Haplotypes are indexed from 1.
#' @param to_recipient last recipient haplotype included if creating a partial forward table.
#'   By default all are included upto the last recipient haplotype.
#'   Haplotypes are indexed from 1.
#'
#' @return
#'   A specialized list of class `kalisForwardTable`.
#'   The elements of the forward table list are:
#'   \describe{
#'     \item{`l`}{denotes the current variant position.
#'       2147483647 indicates a newly created forward table which has not yet been propagated to any variant.}
#'     \item{`alpha`}{is a matrix of rescaled forward probabilities under the Li and Stephens HMM.
#'       Each column of `alpha` corresponds to an independent HMM such that \eqn{\alpha^l_{ji}} is proportional to the probability that haplotype \eqn{j} is copied by haplotype \eqn{i} at variant \eqn{l} and observing haplotype \eqn{i} from variant 1 up through variant \eqn{l}.}
#'     \item{`alpha.f`}{is a vector containing scaling constants needed to continue propagating the HMM (TODO: add ref please see kalis paper for details).}
#'   }
#'
#'   A `kalisForwardTable` also carries with it a checksum key for the parameters `pars` it was provided.
#'   If one attempts to interact a `kalisForwardTable` with a `kalisBackwardTable` with mismatched parameters, an error will be thrown.
#'
#' @seealso
#'   [Forward()] to propagate the newly created `kalisForwardTable`;
#'   [MakeBackwardTable()] to create a corresponding `kalisBackwardTable`.
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
#' fwd
#'
#' # Now ready to run the HMM forward recursions with Forward() ...
#'
#' @aliases print.kalisForwardTable
#' @export
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

  # Create table in C to ensure new memory allocated and no references where not wanted
  .Call(CCall_MakeForwardTable, from_recipient, to_recipient, pars$sha256)
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

  if(x$l == 2147483647L) {
    cat("  Newly created table, currently uninitialised to any variant (ready for Forward function next).\n")
  } else {
    cat(glue("  Current variant = {x$l}"), "\n")
  }
  cat("  Memory consumed: ", ceiling(utils::object.size(x)/1e6)/1e3, "GB.\n")
}



#' Create a `kalisBackwardTable`
#'
#' Allocates the memory for and initializes a backward table.
#'
#' `MakeBackwardTable` initializes a `kalisBackwardTable` object that will store the rescaled backward probabilities for each recipient haplotype.
#' Note: since all other haplotypes loaded by [CacheHaplotypes()] into the kalis cache are taken as potential donor haplotypes, this table will be of size (total haplotypes in cache)x(number of recipients).
#' Also note that this table object is linked to a given set of Li and Stephens hidden Markov model parameters, created by [Parameters()], to ensure consistency for any given HMM run.
#'
#' The returned `kalisBackwardTable` is ready to be propagated to a given target variant with the function [Backward()].
#'
#' Since each column corresponds to an independent Li and Stephens hidden Markov model (ie for each recipient), it is possible to create a partial backward table object which corresponds to a subset of recipients using the `from_recipient` and `to_recipient` arguments.
#'
#' @param pars a `kalisParameters` object specifying the genetics parameters to be associated with this backward table.
#'   These parameters can be set up by using the [Parameters()] function.
#' @param from_recipient first recipient haplotype included if creating a partial backward table.
#'   By default all are included from the first recipient haplotype.
#'   Haplotypes are indexed from 1.
#' @param to_recipient last recipient haplotype included if creating a partial backward table.
#'   By default all are included upto the last recipient haplotype.
#'   Haplotypes are indexed from 1.
#'
#' @return
#'   A specialized list of class \code{kalisBackwardTable}.
#'   The elements of the backward table list are:
#'   \describe{
#'     \item{`l`}{denotes the current variant position.
#'       2147483647 indicates a newly created backward table which has not yet been propagated to any variant.}
#'     \item{`beta`}{is a matrix of rescaled backward probabilities under the Li and Stephens HMM.
#'       Each column of `beta` corresponds to an independent HMM such that \eqn{\beta^l_{ji}} is proportional to the probability of observing haplotype \eqn{i} from variant \eqn{l+1} up through variant \eqn{L} given that haplotype \eqn{j} is copied by haplotype \eqn{i} at variant \eqn{l}.}
#'     \item{`beta.g`}{is a vector containing scaling constants needed to continue propagating the HMM (TODO: add ref please see kalis paper for details).}
#'     \item{`beta.theta`}{boolean indicator for whether the matrix `beta` is currently in so-called "beta-theta space" or not.  See [Backward()] for details about "beta-theta space" which is specified during the run.}
#'   }
#'
#'   A `kalisBackwardTable` also carries with it a checksum key for the parameters `pars` it was provided.
#'   If one attempts to interact a `kalisBackwardTable` with a `kalisForwardTable` with mismatched parameters, an error will be thrown.
#'
#' @seealso
#'   [Backward()] to propagate the newly created `kalisBackwardTable`;
#'   [MakeForwardTable()] to create a corresponding `kalisForwardTable`.
#'
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
#' bck <- MakeBackwardTable(pars)
#' bck
#'
#' # Now ready to run the HMM backward recursions with Backward() ...
#'
#' @aliases print.kalisBackwardTable
#' @export
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

  # Create table in C to ensure new memory allocated and no references where not wanted
  .Call(CCall_MakeBackwardTable, from_recipient, to_recipient, pars$sha256)
}

#' @export
print.kalisBackwardTable <- function(x, ...) {
  if(!("kalisBackwardTable" %in% class(x)))
    stop("Not a kalisBackwardTable object")

  d <- dim(x$beta)

  if(d[1]==d[2]) {
    cat(glue("Full Backward Table object for {d[1]} haplotypes, in {ifelse(x$beta.theta, 'beta-theta', 'rescaled probability')} space."), "\n")
  } else {
    cat(glue("Partial Backward Table object for {d[1]} haplotypes, in {ifelse(x$beta.theta, 'beta-theta', 'rescaled probability')} space."), "\n")
    cat(glue("  Recipients {x$from_recipient} to {x$to_recipient}"), "\n")
  }

  if(x$l == 2147483647L) {
    cat("  Newly created table, currently uninitialised to any variant (ready for Backward function next).\n")
  } else {
    cat(glue("  Current variant = {x$l}"), "\n")
  }
  cat("  Memory consumed: ", ceiling(utils::object.size(x)/1e6)/1e3, "GB.\n")
}

#' Copy Forward/Backward tables
#'
#' Copies the contents of one forward/backward table into another.
#'
#' The core code in kalis operates on forward and backward tables at a very low level, both for speed (using low level CPU vector instructions) but also to avoid unnecessary memory copies since these tables will tend to be very large in serious genetics applications.
#' As a result, if you attempt to copy a table in the standard idomatic way in R:
#' ```
#' fwd2 <- fwd
#' ```
#' then these two variables merely point to the *same* table: running the forward algorithm on `fwd` would result in `fwd2` also changing.
#'
#' This function is therefore designed to enable explicit copying of tables.
#'
#' @param to a `kalisForwardTable` or `kalisBackwardTable` object which is to be copied into.
#' @param from a `kalisForwardTable` or `kalisBackwardTable` object which is to be copied from.
#'
#' @seealso
#'   [MakeForwardTable()], [MakeBackwardTable()] to create tables which can be copied into.
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
#' # Move to variant 10
#' Forward(fwd, pars, 10)
#' fwd
#'
#' # This does **NOT** work as intended:
#' fwd2 <- fwd
#' Forward(fwd, pars, 20)
#'
#' # Both tables are now at variant 20
#' fwd
#' fwd2
#'
#' # Instead, to copy we create another table and use this function
#' fwd2 <- MakeForwardTable(pars)
#' CopyTable(fwd2, fwd)
#'
#' # Now propagate fwd again
#' Forward(fwd, pars, 30)
#'
#' # But notice this fwd2 wasn't affected
#' fwd
#' fwd2
#'
#' @importFrom lobstr ref
#' @export
CopyTable <- function(to, from) {
  if(!("kalisForwardTable" %in% class(to)) && !("kalisBackwardTable" %in% class(to))) {
    stop("The to argument is not a valid forward or backward table.")
  }
  if(!("kalisForwardTable" %in% class(from)) && !("kalisBackwardTable" %in% class(from))) {
    stop("The from argument is not a valid forward or backward table.")
  }

  type <- NULL
  if("kalisForwardTable" %in% class(from))
    type <- "fwd"
  else if("kalisBackwardTable" %in% class(from))
    type <- "bck"

  if(is.null(type)) {
    stop("Error identifying type of from table")
  } else if(type == "fwd" && !("kalisForwardTable" %in% class(to))) {
    stop("type mismatch: from is a forward table, but to is a backward.")
  } else if(type == "bck" && !("kalisBackwardTable" %in% class(to))) {
    stop("type mismatch: from is a backward table, but to is a forward.")
  }

  if(to$from_recipient != from$from_recipient) {
    stop(glue("from table starts with recipient {from$from_recipient}, but to table starts with recipient {to$from_recipient}"))
  }
  if(to$to_recipient != from$to_recipient) {
    stop(glue("from table ends with recipient {from$to_recipient}, but to table ends with recipient {to$to_recipient}"))
  }

  if(to$pars.sha256 != from$pars.sha256) {
    stop("The two tables provided were created with different parameter values (SHA-256 mismatch).")
  }

  if(all(ref(from) == ref(to))) {
    stop("from and to are pointing to the same memory space.  Please make a new table before copying.")
  }

  .Call(CCall_CopyFBTable, to, from)
  invisible(NULL)
}

#' Reset Forward/Backward tables
#'
#' Resets the provided forward/backward table so it can propagate from the start/end of the Li and Stephens hidden Markov model.
#'
#' It is *much* faster to reset a forward/backward table rather than remove and make a new one.
#' This function marks a table as reset so that it will be propagated as if freshly allocated.
#'
#' @param tbl a `kalisForwardTable` or `kalisBackwardTable` object
#'   which is to be reset.
#'
#' @seealso [MakeForwardTable()], [MakeBackwardTable()] to create fresh tables.
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
#' # Move to variant 10
#' Forward(fwd, pars, 10)
#' fwd
#'
#' # This does **NOT** work as intended, cannot reverse from variant 10:
#' try(Forward(fwd, pars, 5))
#'
#' # Could reset the table though:
#' ResetTable(fwd)
#' Forward(fwd, pars, 5)
#' fwd
#'
#' @export
ResetTable <- function(tbl) {
  if(!("kalisForwardTable" %in% class(tbl)) && !("kalisBackwardTable" %in% class(tbl))) {
    stop("The tbl argument is not a valid forward or backward table.")
  }

  type <- NULL
  if("kalisForwardTable" %in% class(tbl))
    type <- "fwd"
  else if("kalisBackwardTable" %in% class(tbl))
    type <- "bck"

  if(is.null(type)) {
    stop("Error identifying type of table")
  }

  if(type == "fwd") {
    .Call(CCall_ResetForwardTable, tbl)
  }
  if(type == "bck") {
    .Call(CCall_ResetBackwardTable, tbl)
  }
}
