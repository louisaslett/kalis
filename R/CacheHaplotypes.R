pkgVars <- new.env(parent = emptyenv())
assign("working.dir", '.', envir = pkgVars)
assign("N", NA, envir = pkgVars) # must be integer
assign("L", NA, envir = pkgVars) # must be integer

# We assume haplotypes are stored in the slowest changing dimension
# per the HDF5 spec definition.  This is "row-wise" in the C standard
# spec, or "col-wise" in the rhdf5 spec.

#' Create memory cache of haplotypes
#'
#' Load haplotypes from hard drive or an R matrix to memory.
#'
#' To achieve higher performance, kalis internally represents haplotypes
#' in an efficient raw binary format in memory.  This function will load
#' haplotypes from a file or from a binary R matrix and convert this into kalis'
#' internal format ready for use by the other functions in this package.  Note
#' that only one set of haplotypes can be cached at a time and calling this
#' function twice overwrites cache of haplotypes created by the first function
#' call.
#'
#' At present, only HDF5 is supported natively and there is a vignette giving a
#' simple script to convert from VCF to HDF5.  If there is sufficient demand
#' other formats may be added in future.
#'
#' **R matrix**
#'
#' If supplying an R matrix, it must consist of only 0's or 1's.  The haplotypes
#' should be stored in columns, with loci in rows.  That is, the dimensions
#' should be:
#'
#' (num rows)x(num cols) = (num loci)x(num haplotypes).
#'
#' It is fine to delete this matrix after calling \code{CacheHaplotypes}.
#'
#' **HDF5 format**
#'
#' For HDF5 files, kalis expects a 2-dimensional object named `haps` at the root
#' level of the HDF5 file.  Haplotypes should be stored in the slowest changing
#' dimension as defined in the HDF5 specification (note that different languages
#' treat this as rows or columns).  If the haplotypes are stored in the other
#' dimension then simply set the argument `transpose = TRUE`.
#' If the user is unsure of the convention of the language they used to create
#' the HDF5 file, then the simplest approach is to simply load the data
#' specifying only the HDF5 file name and then confirm that number of haplotypes
#' and their length have not been exchanged in the diagnostic output which kalis
#' prints.
#'
#' @param haps can be the name of a file from which the haplotypes are to be
#'   read, or can be a binary R matrix.
#' @param loci.idx an optional integer vector of which loci to load from `haps`
#'   into the cache, indexed from 1.
#' @param hap.idx an optional integer vector of which haplotypes to load from
#'   `haps` into the cache, indexed from 1.
#' @param format the file format which `file` is stored in, or `"auto"` to
#'   detect format based on the file extension.  Recognised options are `"hdf5"`
#'   and `"vcf"`
#' @param ... options passed on to the specific format transcoding engine.  The
#'   following additional arguments are supported:
#'
#'   1. HDF5
#'       - `transpose` a logical value indicating whether to switch the
#'         interpretation of the slowest changing dimension (hence switching the
#'         number of haplotypes and the length)
#'       - `haps` a string giving the path to a 2-dimensional object in the
#'         HDF5 file specifying the haplotype matrix.  Defaults to `/haps`
#'       - `hdf5.pkg` a string giving the HDF5 R package to use to load the file
#'         from disk.  The packages `rhdf5` (BioConductor) and `hdf5r` (CRAN)
#'         are both supported.  Default is to use `hdf5r` if both packages are
#'         available, with fallback to `rhdf5`.  This should never need to be
#'         specified unless you have both packages but want to force the use of
#'         the `rhdf5` package.
#'   1. R matrix
#'       - `transpose` a logical value indicating whether to switch the
#'         interpretation of the dimensions (hence switching the number of
#'         haplotypes and the length).  By default (`FALSE`), loci are taken to
#'         be in rows with haplotypes in columns (ie a num loci x num haps
#'         matrix)
#'
#' @return Nothing is returned by the function.
#'   However, a status message is output indicating the dimensions of the loaded
#'   data.
#'   This can be especially useful for HDF5 input when you are uncertain about
#'   the orientation of the `haps` object.
#'   If this message shows the dimensions incorrectly ordered, call the function
#'   again with the extra argument `transpose = TRUE`.
#'
#' @seealso \code{\link{ClearHaplotypeCache}} to remove the haplotypes from
#'   the cache and free the memory.
#'
#' @examples
#' # Examples
#' \dontrun{
#' # Load haplotypes to cache from an HDF5 file on disk
#' CacheHaplotypes("myhaps.h5")
#'
#' # If the diagnostic message printed during the above indicates the numbers
#' # of haplotypes and their length are the wrong way around, reload with
#' # argument to transpose
#' CacheHaplotypes("myhaps.h5", transpose = TRUE)
#'
#' # When correct orientation is known, you can avoid diagnostic messages for
#' # running in script files
#' suppressMessages(CacheHaplotypes("myhaps.h5"))
#'
#' # Alternatively, if you have an exotic file format that can be loaded in to
#' # R by other means, then a binary matrix can be supplied.  This example
#' # randomly simulated a binary matrix to illustrate (in reality load from
#' # disk):
#' n.haps <- 100
#' n.loci <- 20000
#' haps <- matrix(sample(0:1, n.haps*n.loci, replace = TRUE),
#'                nrow = n.loci,
#'                ncol = n.haps)
#' CacheHaplotypes(haps)
#' }
#'
#' @export
CacheHaplotypes <- function(haps, loci.idx = NULL, hap.idx = NULL, format = "auto", ...) {
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
    ext <- stringr::str_to_lower(stringr::str_extract(haps, stringr::regex("\\.[0-9a-z]+$")))

    if(format == "hdf5" ||
       (format == "auto" && (ext == ".h5" ||
                             ext == ".hdf5"))) {
      CacheHaplotypes.hdf5(haps, loci.idx, hap.idx, ...)
      cached <- TRUE
    }
    if(format == "vcf" ||
       (format == "auto" && (ext == ".vcf" ||
                             ext == ".gz" && stringr::str_to_lower(stringr::str_extract(haps, stringr::regex("\\.[0-9a-z]\\.[0-9a-z]+$"))) == ".vcf.gz"))) {
      message("Native import of VCF is currently not supported.  Please see the vignettes for a simple script to convert VCF to HDF5.")
      cached <- TRUE
    }

    if(!cached) {
      stop("Unable to identify file format.")
    }
  } else if(is.matrix(haps) || is.data.frame(haps)) {
    if(!is.integer(haps)) {
      stop("Haplotype matrix must be integer.")
    }
    if(any(!(haps == 1 || haps == 0))) {
      stop("All entries in haplotype matrix must be zero or one.")
    }
    CacheHaplotypes.matrix(haps, loci.idx, hap.idx, ...)
  } else {
    stop("Unable to load haplotypes from the provided object.")
  }
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
#' Retrieve haplotypes from the memory cache, converting the raw binary into
#' a simple integer 0/1 for inspection and use in R.
#'
#' To achieve higher performance, kalis internally represents haplotypes
#' in an efficient raw binary format in memory which cannot be directly viewed
#' or manipulated in R.  This function enables you to copy whole or partial
#' haplotypes out of this low-level format and into a standard R
#' matrix of 0's and 1's.
#'
#' @param hap.idx which haplotypes to retrieve from the cache, specified as a
#'   (vector) index.  This enables specifying haplotypes by offset in the order
#'   they were loaded into the cache (from 1 to the number of haplotypes).
#' @param loci.idx which loci to retrieve from the cache, specified as a
#'   (vector) index.  This enables specifying loci by offset in the order
#'   they were loaded into the cache (from 1 to the number of loci).
#'   Cannot be specified at the same time as `loci.ids`
#'
#' @return A matrix of 0/1 integers with `length` rows and
#'   `length(ids)` columns, such that haplotypes appear in columns.
#'
#' @seealso \code{\link{CacheHaplotypes}} to fill the memory cache with
#'   haplotypes.
#'
#' @examples
#' # Examples
#' \dontrun{
#' # For the purposes of an example, fill the cache with random haplotypes ...
#' n.haps <- 100
#' n.loci <- 20000
#' haps <- matrix(sample(0:1, n.haps*n.loci, replace = TRUE),
#'                nrow = n.loci,
#'                ncol = n.haps)
#' CacheHaplotypes(haps)
#'
#' # ... and confirm we can read a chosen portion back.  Try to read back
#' # the 10th and 11th haplotypes from loci 50 to 250 inclusive
#' res <- QueryCache(10:11, 50, 201)
#' all(res == haps[50:250, 10:11])
#' }
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
      res[i,] <- as.integer(intToBits(QueryCache2_loc(loci.idx[i]-1)))[hap.idx]
    }
  } else {
    for(i in 1:length(hap.idx)) {
      res[,i] <- as.integer(intToBits(QueryCache2_ind(hap.idx[i]-1)))[loci.idx]
    }
  }
  res
}



#' Remove all cached haplotypes and free memory
#'
#' Remove all the haplotypes that were cached by a previous caching call and
#' free up the memory that was allocated for future use.
#'
#' To achieve higher performance, kalis internally represents haplotypes
#' in an efficient raw binary format in memory which cannot be directly viewed
#' or manipulated in R.  In particular, this cache sits outside R's memory
#' management and will never be garbage collected (unless R is quit or the
#' package is unloaded).  Therefore, this function is provided to enable freeing
#' the memory used by this cache.
#'
#' @return Nothing is returned.
#'
#' @seealso \code{\link{CacheHaplotypes}} to create a haplotype cache.
#'
#' @examples
#' # Examples
#' \dontrun{
#' ClearHaplotypeCache()
#' }
#'
#' @export
ClearHaplotypeCache <- function() {
  assign("N", NA, envir = pkgVars)
  assign("L", NA, envir = pkgVars)
  ClearHaplotypeCache2()
}



.onUnload <- function(libpath) {
  assign("N", NA, envir = pkgVars)
  assign("L", NA, envir = pkgVars)
  ClearHaplotypeCache()
  library.dynam.unload("kalis", libpath)
}
