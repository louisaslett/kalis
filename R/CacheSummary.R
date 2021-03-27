#' Retrieve information about the haplotype cache
#'
#' @return
#'   \code{CacheSummary()} prints information about the current state of the kalis
#'   cache.
#'
#'   \code{N()} returns the number of haplotypes currently in the kalis cache, or
#'   \code{NULL} if the cache is empty.
#'
#'   \code{L()} returns the number of variants currently in the kalis cache, or
#'   \code{NULL} if the cache is empty.
#'
#' @export
CacheSummary <- function() {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    cat("Cache currently empty.")
  } else {
    L <- get("L", envir = pkgVars)
    cat(glue("Cache currently loaded with {N} haplotypes, each with {L} variants."), "\n")
    alignment <- 4*(.Call(CCall_VectorBitWidth)/32);
    while(alignment < .Machine$sizeof.pointer) { # POSIX alignment must be at least sizeof(void*)
      alignment <- alignment*2;
    }
    x <- alignment/4
    cat(glue("  Memory consumed ≈ {signif((L*ceiling((N/32.0)/x)*x*4)/1073741824.0, 4)} GB."), "\n")
  }
}

# #' @describeIn NandL Retrieve N
#' @rdname CacheSummary
#' @export
N <- function() {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    return(NULL)
  }
  N
}

# #' @describeIn NandL Retrieve L
#' @rdname CacheSummary
#' @export
L <- function() {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    return(NULL)
  }
  get("L", envir = pkgVars)
}
