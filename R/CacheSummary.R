#' Retrieve information about the haplotype cache
#'
#' @return
#'   `CacheSummary()` prints information about the current state of the kalis cache.
#'     Also invisibly returns a vector giving the dimensions of the cached haplotype data (num variants, num haplotypes), or `NULL` if the cache is empty.
#'
#'   `N()` returns the number of haplotypes currently in the kalis cache, or `NULL` if the cache is empty.
#'
#'   `L()` returns the number of variants currently in the kalis cache, or `NULL` if the cache is empty.
#'
#' @examples
#' # First fill the cache with the toy data included in the package
#' data("SmallHaps")
#' CacheHaplotypes(SmallHaps)
#'
#' # View full summary
#' CacheSummary()
#'
#' # Also note the invisible return
#' dims <- CacheSummary()
#' dims
#'
#' # Get just numbers of haplotypes and variants separately
#' N()
#' L()
#'
#' @importFrom prettyunits pretty_bytes
#' @export
CacheSummary <- function() {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    cat("Cache currently empty.")
    invisible(NULL)
  } else {
    L <- get("L", envir = pkgVars)
    cat(glue("Cache currently loaded with {N} haplotypes, each with {L} variants."), "\n")
    alignment <- 4*(.Call(CCall_VectorBitWidth)/32);
    while(alignment < .Machine$sizeof.pointer) { # POSIX alignment must be at least sizeof(void*)
      alignment <- alignment*2;
    }
    x <- alignment/4
    cat(glue("  Memory consumed: {prettyunits::pretty_bytes(L*ceiling((N/32.0)/x)*x*4)}."), "\n")
    invisible(c(L(), N()))
  }
}

#' @rdname CacheSummary
#' @export
N <- function() {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    return(NULL)
  }
  N
}

#' @rdname CacheSummary
#' @export
L <- function() {
  N <- get("N", envir = pkgVars)
  if(anyNA(N)) {
    return(NULL)
  }
  get("L", envir = pkgVars)
}
