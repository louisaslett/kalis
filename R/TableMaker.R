#' Create a \code{kalisForwardTable}
#'
#' Allocates the memory for and initializes a forward table.
#'
#' \code{MakeForwardTable} initializes a \code{kalisForwardTable} object that will
#' store the rescaled forward probabilities for each recipient haplotype.  Note: since all
#' other haplotypes loaded by \code{\link{CacheHaplotypes}} into the kalis cache
#' are taken as potential donor haplotypes, this
#' table will be of size (total # haplotypes in cache)x(# recipients).
#' Also note that this table object is linked to a given set of Li and Stephens hidden
#' Markov model parameters, created by \code{\link{Parameters}}, to ensure consistency
#' for any given HMM run.
#'
#' The returned \code{kalisForwardTable} is ready
#' to be propagated to a given target variant with the function
#' \code{\link{Forward}}.
#'
#' Note: each column corresponds to an independent Li and Stephens
#' hidden Markov model (ie for each recipient).  Therefore, it is possible to
#' create a partial forward table object which
#' corresponds to a subset of recipients using the \code{from_recipient} and
#' \code{to_recipient} arguments.
#'
#' @param pars a \code{kalisParameters} object specifying the genetics
#'   parameters to be associated with this forward table.  These parameters can
#'   be set up by using the \code{\link{Parameters}} function.
#' @param from_recipient first recipient haplotype included if creating a partial forward
#'   table.  By default includes from the first recipient haplotype.  Haplotypes are
#'   indexed from 1.
#' @param to_recipient last recipient haplotype included if creating a partial forward
#'   table.  By default includes to the last recipient haplotype.  Haplotypes are
#'   indexed from 1.
#'
#' @return
#'   A specialized list of class \code{kalisForwardTable}.
#'   The elements of the forward table list are:
#'   \describe{
#'     \item{\code{l}}{denotes the current variant position.  Zero indicates a
#'       newly created forward table which has not yet been propagated to any variant.}
#'     \item{\code{alpha}}{is a matrix of rescaled forward probabilities under
#'       the Li and Stephens HMM.
#'       Each column of \code{alpha} corresponds to an independent HMM such
#'       that \eqn{\alpha^l_{ji}} is proportional to the probability that
#'       haplotype \eqn{j} is copied by haplotype \eqn{i} at variant \eqn{l} and
#'       observing haplotype \eqn{i} from variant 1 up through variant \eqn{l}.}
#'     \item{\code{alpha.f}}{is a vector containing scaling constants needed to
#'       continue propagating the HMM (TODO: add ref please see kalis paper for details).}
#'   }
#'
#'   A \code{kalisForwardTable} also carries with it a checksum key for the
#'   parameters \code{pars} it was provided.
#'   If one attempts to interact a \code{kalisForwardTable} with a
#'   \code{kalisBackwardTable} with mismatched parameters, an error will be
#'   thrown.
#'
#'
#' @seealso
#'   \code{\link{Forward}} to propagate the newly created
#'   \code{kalisForwardTable}.
#'   \code{\link{MakeBackwardTable}} to create a corresponding
#'   \code{kalisBackwardTable}.
#'
#' @examples
#' # Examples
#' \dontrun{
#' # Load the toy haplotype example and set toy parameters
#' CacheHaplotypes(SmallHaplotypes)
#' pars <- Parameters(rho = CalcRho(cM=SmallMap))
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
#' # which will trigger initialisation and propagation from the first variant
#' # For example, initialise and propagate forward to variant 10:
#' Forward(fwd, pars, 10, nthreads = 8)
#' }
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
    cat("  Newly created table, currently uninitialised to any variant (ready for Forward function next).\n")
  } else {
    cat(glue("  Current variant = {x$l}"), "\n")
  }
  cat("  Memory consumed: ", ceiling(utils::object.size(x)/1e6)/1e3, "GB.\n")
}



#' Create a \code{kalisBackwardTable}
#'
#' Allocates the memory for and initializes a backward table.
#'
#' \code{MakeBackwardTable} initializes a \code{kalisBackwardTable} object that will
#' store the rescaled backward probabilities for each recipient haplotype.  Note: since all
#' other haplotypes loaded by \code{\link{CacheHaplotypes}} into the kalis cache
#' are taken as potential donor haplotypes, this
#' table will be of size (total # haplotypes in cache)x(# recipients).
#' Also note that this table object is linked to a given set of Li and Stephens hidden
#' Markov model parameters, created by \code{\link{Parameters}}, to ensure consistency
#' for any given HMM run.
#'
#' The returned \code{kalisBackwardTable} is ready
#' to be propagated to a given target variant with the function
#' \code{\link{Backward}}.
#'
#' Note: each column corresponds to an independent Li and Stephens
#' hidden Markov model (ie for each recipient).  Therefore, it is possible to
#' create a partial backward table object which
#' corresponds to a subset of recipients using the \code{from_recipient} and
#' \code{to_recipient} arguments.
#'
#' @param pars a \code{kalisParameters} object specifying the genetics
#'   parameters to be associated with this backward table.  These parameters can
#'   be set up by using the \code{\link{Parameters}} function.
#' @param from_recipient first recipient haplotype included if creating a partial backward
#'   table.  By default includes from the first recipient haplotype.  Haplotypes are
#'   indexed from 1.
#' @param to_recipient last recipient haplotype included if creating a partial backward
#'   table.  By default includes to the last recipient haplotype.  Haplotypes are
#'   indexed from 1.
#'
#' @return
#'   A specialized list of class \code{kalisBackwardTable}.
#'   The elements of the backward table list are:
#'   \describe{
#'     \item{\code{l}}{denotes the current variant position.  2147483647 indicates a
#'       newly created backward table which has not yet been propagated to any variant.}
#'     \item{\code{beta}}{is a matrix of rescaled backward probabilities under
#'       the Li and Stephens HMM.
#'       Each column of \code{beta} corresponds to an independent HMM such
#'       that \eqn{\beta^l_{ji}} is proportional to the probability of observing
#'       haplotype \eqn{i} from variant \eqn{l+1} up through variant \eqn{L} given
#'       that haplotype \eqn{j} is copied by haplotype \eqn{i} at variant
#'       \eqn{l}.}
#'     \item{\code{beta.g}}{is a vector containing scaling constants needed to
#'       continue propagating the HMM (TODO: add ref please see kalis paper for details).}
#'     \item{\code{beta.theta}}{boolean indicator for whether the matrix beta
#'       is currently in so-called beta-theta space or not.  See \code{\link{Backward}}
#'       for details about beta-theta space which is specified during the run.}
#'   }
#'
#'   A \code{kalisBackwardTable} also carries with it a checksum key for the
#'   parameters \code{pars} it was provided.
#'   If one attempts to interact a \code{kalisBackwardTable} with a
#'   \code{kalisForwardTable} with mismatched parameters, an error will be
#'   thrown.
#'
#' @seealso \code{\link{Backward}} to propagate the newly created \code{kalisBackwardTable}.
#'
#' @examples
#' \dontrun{
#' # Load the toy haplotype example and set toy parameters
#' CacheHaplotypes(SmallHaplotypes)
#' pars <- Parameters(rho = CalcRho(cM=SmallMap))
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
#' # For example, initialise and propagate backward to variant 10:
#' Backward(bck, pars, 10, nthreads = 8)
#' }
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
  delN <- to_recipient-from_recipient+1

  bck <- list()

  # Define core table, duplicating where relevant to ensure unique to this forward table
  bck$beta           <- duplicate(matrix(0, N, delN))
  bck$beta.g         <- duplicate(rep(0, delN))
  bck$beta.theta     <- duplicate(FALSE)
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
    cat(glue("Full Backward Table object for {d[1]} haplotypes, in {ifelse(x$beta.theta, 'beta-theta', 'rescaled probability')} space."), "\n")
  } else {
    cat(glue("Partial Backward Table object for {d[1]} haplotypes, in {ifelse(x$beta.theta, 'beta-theta', 'rescaled probability')} space."), "\n")
    cat(glue("  Recipients {x$from_recipient} to {x$to_recipient}"), "\n")
  }

  if(x$l == 2147483647) {
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
#' The core code in kalis operates on forward and backward tables at a very low
#' level, both for speed (using low level CPU vector instructions) but also to
#' avoid unnecessary memory copies since these tables will tend to be very large
#' in serious genetics applications.  As a result, if you attempt to copy a table
#' in the standard idomatic way in R:
#'
#' \code{fwd2 <- fwd}
#'
#' then these two variables merely point to the *same* table: running the
#' forward algorithm on \code{fwd} would result in \code{fwd2} also changing.
#'
#' This function is therefore designed to enable explicit copying of tables.
#'
#' @param to a \code{kalisForwardTable} or \code{kalisBackwardTable} object
#'   which is to be copied into.
#' @param from a \code{kalisForwardTable} or \code{kalisBackwardTable} object
#'   which is to be copied from.
#'
#' @seealso \code{\link{MakeForwardTable}}, \code{\link{MakeForwardTable}} to
#' create tables which can be copied into.
#'
#' @examples
#' # Examples
#' \dontrun{
#' # Create a forward table for the hidden Markov model incorporating all
#' # recipient and donor haplotypes
#' fwd <- MakeForwardTable()
#'
#' # Propagate forward to variant 10:
#' Forward(fwd, pars, 10, nthreads = 8)
#'
#' # This does **NOT** work as intended:
#' fwd2 <- fwd
#' Forward(fwd, pars, 20, nthreads = 8)
#'
#' # Both tables are now at variant 20
#' fwd
#' fwd2
#'
#' # Instead, to copy we create another table and use this function
#' fwd2 <- MakeForwardTable()
#' CopyTable(fwd2, fwd)
#' }
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

  if(type == "fwd") {
    CopyForwardTable_cpp(to, from)
  }
  if(type == "bck") {
    CopyBackwardTable_cpp(to, from)
  }
}
