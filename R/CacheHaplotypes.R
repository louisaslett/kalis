pkgVars <- new.env(parent = emptyenv())
assign("working.dir", '.', envir = pkgVars)
assign("N", NA, envir = pkgVars) # must be integer
assign("L", NA, envir = pkgVars) # must be integer

# We assume haplotypes are stored in the slowest changing dimension
# per the HDF5 spec definition.  This is "row-wise" in the C standard
# spec, or "col-wise" in the rhdf5 spec.

#' Load haplotypes into optimised package cache
#'
#' Load haplotypes from hard drive or an R matrix into an optimised kalis package memory cache (overwrites any previous load).
#'
#' To achieve higher performance, kalis internally represents haplotypes in an efficient raw binary format in memory.
#' This function will load haplotypes from a file or from a binary R matrix and convert this into kalis' internal format ready for use by the other functions in this package.
#' Note that only one set of haplotypes can be cached at a time and calling this function twice overwrites cache of haplotypes created by the first function call.
#'
#' Including singletons (variants where there is only one 1 or only one 0) in the loaded haplotypes can lead to numerical instability and columns of `NaN`s in the resulting forward and backward tables when `mu` (see [Parameters()]) is small.
#' Thus, kalis throws a warning when loaded haplotypes contain singletons.
#'
#' At present, hap.gz and hdf5 are supported natively, see the Examples section below showing for how to convert from a VCF/BCF to hap.gz with one `bcftools` command.
#'
#'
#' **hap.gz format**
#'
#' This is the HAP/LEGEND/SAMPLE format used by IMPUTE2 and SHAPEIT.
#' Only the `.hap.gz` file is required for loading with `CacheHaplotypes`, though the `.legend.gz` file can speed up reading the haplotypes.
#' See \url{http://samtools.github.io/bcftools/bcftools.html#convert} for more details on this format.
#'
#'
#' **R matrix**
#'
#' If supplying an R matrix, it must consist of only 0's or 1's.
#' The haplotypes should be stored in columns, with variants in rows.
#' That is, the dimensions should be:
#'
#' (num rows)x(num cols) = (num variants)x(num haplotypes).
#'
#' It is fine to delete this matrix from R after calling [CacheHaplotypes()].
#'
#'
#' **HDF5 format**
#'
#' For HDF5 files, kalis expects a 2-dimensional object named `haps` at the root level of the HDF5 file.
#' Haplotypes should be stored in the slowest changing dimension as defined in the HDF5 specification (note that different languages treat this as rows or columns).
#' If the haplotypes are stored in the other dimension then simply set the argument `transpose = TRUE`.
#' If the user is unsure of the convention of the language they used to create the HDF5 file, then the simplest approach is to simply load the data specifying only the HDF5 file name and then confirm that number of haplotypes and their length have not been exchanged in the diagnostic output which kalis prints.
#'
#'
#'
#' @references
#' Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation of the Li & Stephens model for local ancestry inference in R", *BMC Bioinformatics*, **25**(1). Available at: \doi{10.1186/s12859-024-05688-8}.
#'
#' @param haps can be the name of a file from which the haplotypes are to be read, or can be an R matrix containing only 0/1s.
#'   See Details section for supported file types.
#' @param loci.idx an optional vector of indices specifying the variants to load into the cache, indexed from 1.
#' @param hap.idx an optional vector of indices specifying the haplotypes to load into the cache, indexed from 1.
#' @param warn.singletons a logical, if `FALSE`, suppress warning that singletons (variants where there is only one 1 or only one 0) are present in the loaded `haps`.
#' @param format the file format that `haps` is stored in, or `"auto"` to detect the format based on the file extension.
#'   Recognised options are `"hapgz"` (format used by IMPUTE2 and SHAPEIT) or `"hdf5"` (custom).
#'   See Details section for more information, and for easy conversion from VCF/BCF and other formats see the Examples section.
#' @param ... format specific options for reading in `haps`.
#'   Supported optional arguments for each format are:
#'
#'   1. For `"hapgz"`
#'       - `legendgz.file` a string for faster loading: a `.legend.gz` file can be supplied and will be used to more efficiently determine the number of variants in the `.hap.gz` file
#'       - `L` an integer for faster loading: the number of variants in the `.hap.gz` file can be directly provided
#'       - `N` an integer for faster loading: the number of haplotypes in the `.hap.gz` file can be directly provided
#'   1. For `"hdf5"`
#'       - `transpose` a logical, if `TRUE`, switch the interpretation of rows and columns in `haps`: hence switching the number of haplotypes and the number of variants (the HDF5 specification does not prescribe row/column interpretation, only defining the slowest changing dimension as 'first').
#'         Defaults to `FALSE`.
#'       - `haps.path` a string giving the path to a 2-dimensional object in the HDF5 file specifying the haplotype matrix.
#'         Defaults to `/haps`
#'       - `hdf5.pkg` a string giving the HDF5 R package to use to load the file from disk.
#'         The packages `rhdf5` (BioConductor) and `hdf5r` (CRAN) are both supported.
#'         Default is to use `hdf5r` if both packages are available, with fallback to `rhdf5`.
#'         This should never need to be specified unless you have both packages but want to force the use of the `rhdf5` package.
#'   1. R matrix
#'       - `transpose` a logical, if `TRUE`, switch the interpretation of rows and columns in `haps`: hence switching the number of haplotypes and the number of variants.
#'         Defaults to `FALSE`, meaning variants are taken to be in rows with haplotypes in columns (ie a num variants x num haplotypes matrix)
#'
#' @return A vector giving the dimensions of the cached haplotype data is invisibly returned (num variants, num haplotypes).
#'   It is highly recommended that you run [CacheSummary()] after `CacheHaplotypes`, especially if you are uncertain about the interpretation of rows and columns in `haps`.
#'   If [CacheSummary()] shows that the number of haplotypes and variants are reversed, try calling `CacheHaplotypes` again with the extra argument `transpose = TRUE`.
#'
#' @seealso [CacheSummary()] for a list detailing the current cache status;
#'   [QueryCache()] to copy the haplotypes in the `kalis` cache into an R matrix;
#'   [ClearHaplotypeCache()] to remove the haplotypes from the cache and free the memory;
#'   [N()] for the number of haplotypes cached;
#'   [L()] for the number of variants cached.
#'
#' @examples
#' \dontrun{
#' # If starting from a VCF/BCF first use bcftools to convert to
#' # HAP/SAMPLE/LEGEND format (bcftools can take in several starting formats)
#' # See http://samtools.github.io/bcftools/bcftools.html#convert
#' system("bcftools convert -h my.vcf.gz")
#' CacheHaplotypes("my.hap.gz")
#' CacheSummary()
#' }
#'
#' # If starting directly from a hap.gz file on disk (HAP/LEGEND/SAMPLE format)
#' \dontrun{
#' CacheHaplotypes("my.hap.gz")
#' }
#' # For example, to load the mini example built into the package:
#' CacheHaplotypes(system.file("small_example/small.hap.gz", package = "kalis"))
#' CacheSummary()
#'
#'
#' # If starting from an HDF5 file on disk
#' \dontrun{
#' CacheHaplotypes("my.h5")
#' }
#' # For example, to load the mini example built into the package:
#' CacheHaplotypes(system.file("small_example/small.h5", package = "kalis"))
#' CacheSummary()
#'
#'
#' # If CacheSummary() indicates that the numbers of haplotypes and variants are
#' # the wrong way around, reload with argument transpose set to TRUE
#' \dontrun{
#' CacheHaplotypes("myhaps.h5", transpose = TRUE)
#' CacheSummary()
#' }
#'
#'
#' # Alternatively, if you have an exotic file format that can be loaded in to R
#' # by other means, then a binary matrix can be supplied.  This example
#' # randomly simulates a binary matrix to illustrate.
#' n.haps <- 100
#' n.vars <- 200
#' haps <- matrix(sample(0:1, n.haps*n.vars, replace = TRUE),
#'                nrow = n.vars, ncol = n.haps)
#' CacheHaplotypes(haps)
#' # For example, to load the mini example built into the package:
#' data("SmallHaps")
#' CacheHaplotypes(SmallHaps)
#'
#' @export
CacheHaplotypes <- function(haps, loci.idx = NULL, hap.idx = NULL, warn.singletons = TRUE, format = "auto", ...) {
  if(!identical(loci.idx, NULL)) {
    if(!testAtomicVector(loci.idx, any.missing = FALSE, unique = TRUE)) {
      CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkAtomicVector(loci.idx, any.missing = FALSE, unique = TRUE)}"))
    }
    if(!testIntegerish(loci.idx, lower = 1, sorted = TRUE)) {
      CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, sorted = TRUE)}"))
    }
  }

  if(!identical(hap.idx, NULL)) {
    if(!testAtomicVector(hap.idx, any.missing = FALSE, unique = TRUE)) {
      CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkAtomicVector(hap.idx, any.missing = FALSE, unique = TRUE)}"))
    }
    if(!testIntegerish(hap.idx, lower = 1, sorted = TRUE)) {
      CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, sorted = TRUE)}"))
    }
  }

  if(testString(haps)) {
    if(!testFile(haps, access = "r")) {
      CacheHaplotypes.err(checkFile(haps, access = "r"))
    }

    cached <- FALSE
    ext.pos <- regexpr("\\.[[:alnum:]]+$", haps)
    ext <- if(ext.pos > -1L) { regmatches(haps, ext.pos) } else { "" }
    ext2.pos <- regexpr("\\.[[:alnum:]]+\\.[[:alnum:]]+$", haps)
    ext2 <- if(ext2.pos > -1L) { regmatches(haps, ext2.pos) } else { "" }

    if(format == "hdf5" ||
       (format == "auto" && (ext == ".h5" ||
                             ext == ".hdf5"))) {
      CacheHaplotypes.hdf5(haps, loci.idx, hap.idx, ...)
      cached <- TRUE
    }
    if(format == "vcf" ||
       (format == "auto" && (ext == ".vcf" ||
                             ext == ".gz" && !is.na(ext2) && ext2 == ".vcf.gz"))) {
      message("Native import of VCF is currently not supported.  Please see the vignettes for a simple script to convert VCF to a .hap.gz or HDF5 file.")
      cached <- TRUE
    }
    if(format == "hapgz" ||
       (format == "auto" && (ext == ".gz" && !is.na(ext2) && ext2 == ".hap.gz"))) {
      CacheHaplotypes.hapgz(haps, loci.idx, hap.idx, ...)
      cached <- TRUE
    }

    if(!cached) {
      stop("Unable to identify file format.")
    }
  } else if(is.matrix(haps) || is.data.frame(haps)) {
    if(!is.integer(haps)) {
      stop("Haplotype matrix must be integer.")
    }
    if(any(!(haps == 1 | haps == 0))) {
      stop("All entries in haplotype matrix must be zero or one.")
    }
    CacheHaplotypes.matrix(haps, loci.idx, hap.idx, ...)
  } else {
    stop("Unable to load haplotypes from the provided object.")
  }
  invisible(c(L(), N()))
}



CacheHaplotypes.err <- function(err) {
  if(!is.na(get("N", envir = pkgVars))) {
    stop(glue("{err}  Keeping existing cache."))
  } else {
    stop(err)
  }
}



#' Retrieve haplotypes from memory cache
#'
#' Retrieve haplotypes from the memory cache, converting the raw binary into an integer 0/1 matrix for inspection and use in R.
#'
#' To achieve higher performance, kalis internally represents haplotypes in an efficient raw binary format in memory which cannot be directly viewed or manipulated in R.
#' This function enables you to copy whole or partial views of haplotypes/variants out of this low-level format and into a standard R matrix of 0's and 1's.
#'
#' @references
#' Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation of the Li & Stephens model for local ancestry inference in R", *BMC Bioinformatics*, **25**(1). Available at: \doi{10.1186/s12859-024-05688-8}.
#'
#' @param loci.idx which variants to retrieve from the cache, specified as a (vector) index.
#'   This enables specifying variants by offset in the order they were loaded into the cache (from 1 to the number of variants).
#' @param hap.idx which haplotypes to retrieve from the cache, specified as a (vector) index.
#'   This enables specifying haplotypes by offset in the order they were loaded into the cache (from 1 to the number of haplotypes).
#'
#' @return A matrix of 0/1 integers with `length(loci.idx)` rows and `length(hap.idx)` columns, such that haplotypes appear in columns.
#'
#' @seealso [CacheHaplotypes()] to fill the memory cache with haplotype data.
#'
#' @examples
#' # For the purposes of an example, fill the cache with random haplotypes ...
#' n.haps <- 100
#' n.vars <- 200
#' haps <- matrix(sample(0:1, n.haps*n.vars, replace = TRUE),
#'                nrow = n.vars, ncol = n.haps)
#' CacheHaplotypes(haps)
#'
#' # ... and confirm we can read a chosen portion back.  Try to read back the
#' # 10th and 11th haplotypes from variants 50 to 150 inclusive
#' res <- QueryCache(50:150, 10:11)
#' all(res == haps[50:150, 10:11])
#'
#' @export QueryCache
QueryCache <- function(loci.idx = NULL, hap.idx = NULL) {
  N <- get("N", envir = pkgVars)
  L <- get("L", envir = pkgVars)

  if(!identical(loci.idx, NULL)) {
    if(!testAtomicVector(loci.idx, any.missing = FALSE, unique = TRUE)) {
      CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkAtomicVector(loci.idx, any.missing = FALSE, unique = TRUE)}"))
    }
    if(!testIntegerish(loci.idx, lower = 1, upper = L)) {
      CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, upper = L)}"))
    }
  } else {
    loci.idx <- 1:L
  }

  if(!identical(hap.idx, NULL)) {
    if(!testAtomicVector(hap.idx, any.missing = FALSE, unique = TRUE)) {
      CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkAtomicVector(hap.idx, any.missing = FALSE, unique = TRUE)}"))
    }
    if(!testIntegerish(hap.idx, lower = 1, upper = N)) {
      CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, upper = N)}"))
    }
  } else {
    hap.idx <- 1:N
  }

  # Finally, grab them!
  res <- matrix(integer(1),
                nrow = length(loci.idx),
                ncol = length(hap.idx))

  if(length(loci.idx) < length(hap.idx)) {
    for(i in 1:length(loci.idx)) {
      res[i,] <- as.integer(intToBits(.Call(CCall_QueryCache2_loc, loci.idx[i]-1)))[hap.idx]
    }
  } else {
    for(i in 1:length(hap.idx)) {
      res[,i] <- as.integer(intToBits(.Call(CCall_QueryCache2_ind, hap.idx[i]-1)))[loci.idx]
    }
  }
  res
}



#' Remove all cached haplotypes and free memory
#'
#' Remove all the haplotypes that were cached by a previous caching call and free up the memory that was allocated for future use.
#'
#' To achieve higher performance, kalis internally represents haplotypes in an efficient raw binary format in memory which cannot be directly viewed or manipulated in R, though you can extract a view from it using [QueryCache()].
#' In particular, this cache sits outside R's memory management and will never be garbage collected (unless R is quit or the package is unloaded).
#' Therefore, this function is provided to enable freeing the memory used by this cache.
#'
#' @references
#' Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation of the Li & Stephens model for local ancestry inference in R", *BMC Bioinformatics*, **25**(1). Available at: \doi{10.1186/s12859-024-05688-8}.
#'
#' @return Nothing is returned.
#'
#' @seealso [CacheHaplotypes()] to create a haplotype cache;
#'   [QueryCache()] to view the cache contents in an R matrix.
#'
#' @examples
#' # First fill the cache with the toy data included in the package
#' data("SmallHaps")
#' CacheHaplotypes(SmallHaps)
#'
#' # Verify it is there
#' CacheSummary()
#'
#' # Now clear
#' ClearHaplotypeCache()
#'
#' # Verify it is gone
#' CacheSummary()
#'
#' @export
ClearHaplotypeCache <- function() {
  assign("N", NA, envir = pkgVars)
  assign("L", NA, envir = pkgVars)
  .Call(CCall_ClearHaplotypeCache2)
  invisible(NULL)
}



.onUnload <- function(libpath) {
  assign("N", NA, envir = pkgVars)
  assign("L", NA, envir = pkgVars)
  ClearHaplotypeCache()
  library.dynam.unload("kalis", libpath)
}
