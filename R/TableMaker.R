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
MakeForwardTable <- function(morgan.dist = 0, Ne = 1, gamma = 1, mu = 1e-8, Pi = NULL, from_recipient = 1, to_recipient = Inf) {
  haps <- get("haps", envir = pkgCache)
  if(anyNA(haps)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  N <- length(haps)
  L <- get("hap_size", envir = pkgCache)


  if(!is.numeric(morgan.dist)) {
    stop("morgan.dist must be numeric vector type.")
  }

  if(length(morgan.dist)==1){
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
  if(!is.vector(gamma) || !is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
    stop("gamma must be a positive scalar.")
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

  if(is.null(Pi)){Pi <- 1/(N-1)}

  if(is.data.frame(Pi)) {
    stop("Pi must be a matrix or left as NULL to have uniform copying probabilities of 1/(N-1) for a problem with N recipients, not a data frame.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.vector(Pi) && is.numeric(Pi) && length(Pi) == 1 && Pi == 1/(N-1))) {
    stop("Pi must be a matrix or left as NULL to have uniform copying probabilities of 1/(N-1) for a problem with N recipients.")
  }

  # Compute rho ... this is all we actually want to carry round with the forward
  rho <- c(1-exp(-Ne*morgan.dist^gamma), 1)
  rho <- ifelse(rho<1e-16, 1e-16, rho)


  if(from_recipient>to_recipient) {
    stop("from_recipient must be smaller than to_recipient.")
  }
  if(from_recipient < 1) {
    from_recipient <- 1
  }
  if(to_recipient > N) {
    to_recipient <- N
  }
  delN <- to_recipient-from_recipient+1

  # Define core table, duplicating to ensure unique to this forward table
  fwd <- duplicate(list(alpha          = matrix(0, N, delN),
                        alpha.f        = rep(0, delN),
                        alpha.f2       = rep(0, delN),
                        l              = c(0),
                        from_recipient = from_recipient,
                        to_recipient   = to_recipient))
  # Setup class (could invoke object duplication even withou dupliate above?)
  class(fwd) <- "kalisForwardTable"

  # Add objects we want to explicitly be reference counted last
  fwd$rho <- rho
  fwd$mu <- mu
  fwd$Pi <- Pi

  fwd
}

print.kalisForwardTable <- function(x, ...) {
  if(class(x)!="kalisForwardTable")
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
MakeBackwardTable <- function(morgan.dist = 0, Ne = 1, gamma = 1, mu = 1e-8, Pi = NULL, from_recipient = 1, to_recipient = Inf) {
  haps <- get("haps", envir = pkgCache)
  if(anyNA(haps)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  N <- length(haps)
  L <- get("hap_size", envir = pkgCache)

  if(!is.numeric(morgan.dist)) {
    stop("morgan.dist must be numeric vector type.")
  }

  if(length(morgan.dist)==1){
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
  if(!is.vector(gamma) || !is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
    stop("gamma must be a positive scalar.")
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

  if(is.null(Pi)){Pi <- 1/(N-1)}

  if(is.data.frame(Pi)) {
    stop("Pi must be a matrix or left as NULL to have uniform copying probabilities of 1/(N-1) for a problem with N recipients, not a data frame.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.vector(Pi) && is.numeric(Pi) && length(Pi) == 1 && Pi == 1/(N-1))) {
    stop("Pi must be a matrix or left as NULL to have uniform copying probabilities of 1/(N-1) for a problem with N recipients.")
  }

  # Compute rho ... this is all we actually want to carry round with the forward
  rho <- c(1-exp(-Ne*morgan.dist^gamma), 1)
  rho <- ifelse(rho<1e-16, 1e-16, rho)

  if(from_recipient>to_recipient) {
    stop("from_recipient must be smaller than to_recipient.")
  }
  if(from_recipient < 1) {
    from_recipient <- 1
  }
  if(to_recipient > N) {
    to_recipient <- N
  }
  delN <- to_recipient-from_recipient+1

  bck <- duplicate(list(beta           = matrix(0, N, delN),
                        beta.g         = rep(0, delN),
                        beta.g2        = rep(0, delN),
                        l              = c(2147483647),
                        from_recipient = from_recipient,
                        to_recipient   = to_recipient))

  # Setup class (could invoke object duplication even withou dupliate above?)
  class(bck) <- "kalisBackwardTable"

  # Add objects we want to explicitly be reference counted last
  bck$rho <- rho
  bck$mu <- mu
  bck$Pi <- Pi

  bck
}

print.kalisBackwardTable <- function(x, ...) {
  if(class(x)!="kalisBackwardTable")
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

