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
#' @seealso \code{\link{Forward}} to propagate the newly created table forward
#'   through the genome.
#'
#' @examples
#' # Examples
MakeForwardTable <- function(from_recipient = 1, to_recipient = Inf) {
  seqs <- get("seqs", envir = pkgCache)
  if(anyNA(seqs)) {
    stop("No sequences cached ... cannot determine table size until cache is loaded with CacheAllSequences().")
  }
  N <- length(seqs)
  L <- get("seq_size", envir = pkgCache)

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

  fwd <- duplicate(list(alpha          = matrix(0, N, delN),
                        alpha.f        = rep(0, delN),
                        alpha.f2       = rep(0, delN),
                        l              = c(0),
                        from_recipient = from_recipient,
                        to_recipient   = to_recipient))

  class(fwd) <- "kalisForwardTable"
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
MakeBackwardTable <- function(from_recipient = 1, to_recipient = Inf) {
  seqs <- get("seqs", envir = pkgCache)
  if(anyNA(seqs)) {
    stop("No sequences cached ... cannot determine table size until cache is loaded with CacheAllSequences().")
  }
  N <- length(seqs)
  L <- get("seq_size", envir = pkgCache)

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

  class(bck) <- "kalisBackwardTable"
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

