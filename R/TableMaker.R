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

  list(alpha          = matrix(0, N, delN),
       alpha.f        = rep(0, delN),
       alpha.f2       = rep(0, delN),
       l              = c(0),
       from_recipient = from_recipient,
       to_recipient   = to_recipient)
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

  list(beta           = matrix(0, N, delN),
       beta.g         = rep(0, delN),
       beta.g2        = rep(0, delN),
       l              = c(2147483647),
       from_recipient = from_recipient,
       to_recipient   = to_recipient)
}
