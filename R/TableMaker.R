#' Create a Forward Table
#'
#' Allocates the memory for and sets up a forward table.
#'
#' For performance and numerical stability reasons kalis operates with lag
#' scaled probabilities and these additional scaling factors must be tracked to
#' enable recovery of standard forward probabilities.
#' This utility function allocates memory for the forward table and also sets up
#' the necessary internal tracking of correct scaling factors.
#' However, note that this utility is purely for creating the correct object
#' to track forward computations and upon creation the forward table itself is
#' uninitialised.
#'
#' Therefore, the standard workflow with kalis is to create a forward table
#' using this utility function before passing that to the \code{\link{Forward}}
#' function.  \code{\link{Forward}} will identify uninitialised tables and propagate
#' them correctly from the first locus.
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
#' @return an object of class \code{kalisForwardTable} containing
#'
#' @seealso \code{\link{Forward}} to propagate the newly created table forward
#'   through the genome.
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
#' @export MakeForwardTable
MakeForwardTable <- function(pars, from_recipient = 1, to_recipient = Inf) {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

  if(!("kalisParameters" %in% class(pars))) {
    if("rstudioapi" %in% installed.packages()[, "Package"]) {
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
  fwd$alpha.f2       <- duplicate(rep(0, delN))
  fwd$l              <- duplicate(c(0L))
  fwd$from_recipient <- duplicate(from_recipient)
  fwd$to_recipient   <- duplicate(to_recipient)
  fwd$pars.sha256    <- duplicate(pars$sha256)

  class(fwd) <- c("kalisForwardTable", class(fwd))

  fwd
}

#' @export print.kalisForwardTable
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
  cat("  Memory consumed ≈", ceiling(object.size(x)/1e6)/1e3, "GB.\n")
}

#' Title
#'
#' Short description
#'
#' Detailed description
#'
#' @param from_recipient ...
#' @param to_recipient ...
#'
#' @return Return value
#'
#' @seealso \code{\link{Backward}} to propagate the newly created table forward
#'   through the genome.
#'
#' @examples
#' # Examples
#' @export MakeBackwardTable
MakeBackwardTable <- function(pars, from_recipient = 1, to_recipient = Inf) {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  L <- get("L", envir = pkgVars)

  if(!("kalisParameters" %in% class(pars))) {
    if("rstudioapi" %in% installed.packages()[, "Package"]) {
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

#' @export print.kalisBackwardTable
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
  cat("  Memory consumed ≈", ceiling(object.size(x)/1e6)/1e3, "GB.\n")
}
