CacheHaplotypes.matrix <- function(x, loci.idx, hap.idx, transpose = FALSE) {
  if(!identical(loci.idx, NULL)) {
    if(!transpose) {
      if(!testIntegerish(loci.idx, lower = 1, upper = nrow(x), sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, upper = nrow(x), sorted = TRUE)}"))
      }
    } else {
      if(!testIntegerish(loci.idx, lower = 1, upper = ncol(x), sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, upper = ncol(x), sorted = TRUE)}"))
      }
    }
    loci.idx2 <- loci.idx
  } else {
    if(!transpose) {
      loci.idx2 <- 1:nrow(x)
    } else {
      loci.idx2 <- 1:ncol(x)
    }
  }

  if(!identical(hap.idx, NULL)) {
    if(!transpose) {
      if(!testIntegerish(hap.idx, lower = 1, upper = ncol(x), sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, upper = ncol(x), sorted = TRUE)}"))
      }
    } else {
      if(!testIntegerish(hap.idx, lower = 1, upper = nrow(x), sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, upper = nrow(x), sorted = TRUE)}"))
      }
    }
    hap.idx2 <- hap.idx
  } else {
    if(!transpose) {
      hap.idx2 <- 1:ncol(x)
    } else {
      hap.idx2 <- 1:nrow(x)
    }
  }

  if(identical(hap.idx, NULL) && identical(loci.idx, NULL)) {
    # Do nothing ... no need to subset x
  } else {
    if(!transpose) {
      x <- x[loci.idx2, hap.idx2]
    } else {
      x <- x[hap.idx2, loci.idx2]
    }
  }

  if(!is.na(get("N", envir = pkgVars))) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }

  # Get dimensions of haplotype data
  if(!transpose) {
    N <- dim(x)[2]
    L <- dim(x)[1]
  } else {
    N <- dim(x)[1]
    L <- dim(x)[2]
  }
  assign("N", as.integer(N), envir = pkgVars)
  assign("hap.ids", as.integer(1:N), envir = pkgVars)
  assign("loci.ids", as.integer(1:L), envir = pkgVars)

  # Cache it!
  assign("L", as.integer(CacheHaplotypes_matrix_2(x, N, L, transpose)), envir = pkgVars)
}
