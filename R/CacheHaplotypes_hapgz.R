CacheHaplotypes.hapgz <- function(hapgz.file,
                                  loci.idx,
                                  hap.idx,
                                  legendgz.file = NULL,
                                  L = NULL,
                                  N = NULL) {
  if(!testFile(hapgz.file, access = "r")) {
    CacheHaplotypes.err(checkFile(hapgz.file, access = "r"))
  }

  ext.pos <- regexpr("\\.[[:alnum:]]+\\.[[:alnum:]]+$", hapgz.file)
  ext <- if(ext.pos > -1L) { regmatches(hapgz.file, ext.pos) } else { "" }

  # Figure out L, favouring fastest available option
  if(testCount(L, positive = TRUE)) {
    # use L directly
  } else if(!identical(legendgz.file, NULL)) {
    if(!testFile(legendgz.file, access = "r")) {
      CacheHaplotypes.err(checkFile(legendgz.file, access = "r"))
    }

    L <- .Call(CCall_CacheHaplotypes_hapgz_nlines, legendgz.file) - 1
  } else if(ext == ".hap.gz") {
    legendgz.file <- hapgz.file
    regmatches(legendgz.file, ext.pos) <- ".legend.gz"
    if(!testFile(legendgz.file, access = "r")) {
      CacheHaplotypes.err(glue("Problem with file {legendgz.file} ... {checkFile(legendgz.file, access = 'r')}\n  TIP: specify the location of the legend gzip file using 'legendgz.file' argument, or specify number of variants directly with 'L' argument."))
    }

    L <- .Call(CCall_CacheHaplotypes_hapgz_nlines, legendgz.file) - 1
  } else {
    warning("legendgz.file not specified, cannot auto-detect it, and L not separately specified.  Performing slower scan of haplotype gz file to determine number of loci.\n  TIP: specify the location of the legend gzip file using 'legendgz.file' argument, or specify number of variants directly with 'L' argument to avoid this.")

    L <- .Call(CCall_CacheHaplotypes_hapgz_nlines, hapgz.file)
  }

  # Figure out N, favouring fastest available option
  if(testCount(N, positive = TRUE)) {
    # use N directly
  } else {
    # Figure out N by reading first line of hap.gz
    N <- (.Call(CCall_CacheHaplotypes_hapgz_ncols, hapgz.file) + 1)/2
    if(!testIntegerish(N)) {
      CacheHaplotypes.err("Line length in .hap.gz file is not consistent with haplotype data.")
    }
  }

  # Verify indices
  if(!identical(loci.idx, NULL)) {
    if(!testIntegerish(loci.idx, lower = 1, upper = L, sorted = TRUE)) {
      CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, upper = L, sorted = TRUE)}"))
    }
  } else {
    loci.idx <- 1:L
  }

  if(!identical(hap.idx, NULL)) {
    if(!testIntegerish(hap.idx, lower = 1, upper = N, sorted = TRUE)) {
      CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, upper = N, sorted = TRUE)}"))
    }
  } else {
    hap.idx <- 1:N
  }

  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("N", envir = pkgVars))) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }
  assign("N", as.integer(length(hap.idx)), envir = pkgVars)
  assign("L", as.integer(length(loci.idx)), envir = pkgVars)

  # Cache it!
  .Call(CCall_CacheHaplotypes_hapgz_2, hapgz.file, loci.idx, hap.idx, L, N)

  invisible(NULL)
}


