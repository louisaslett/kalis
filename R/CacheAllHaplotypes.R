pkgVars <- new.env(parent = emptyenv())
assign("working.dir", '.', envir = pkgVars)
assign("N", NA, envir = pkgVars)
assign("L", NA, envir = pkgVars)

# We assume haplotypes are stored in the slowest changing dimension
# per the HDF5 spec definition.  This is "row-wise" in the C standard
# spec, or "col-wise" in the rhdf5 spec.

#' Create memory cache of haplotypes
#'
#' Load haplotypes from hard drive to memory.
#'
#' To achieve higher performance, kalis internally represents haplotypes
#' in an efficient raw binary format in memory.  This function will load haplotypes
#' from a file and convert this into kalis' internal format ready for use
#' by the other functions in this package.
#'
#' At present, only HDF5 is supported natively and there is a vignette giving a simple
#' script to convert from VCF to HDF5.  If there is sufficient demand other formats
#' may be added in future
#'
#' **HDF5 format**
#'
#' For HDF5 files, kalis expects a 2-dimensional object named `haps` at the root
#' level of the HDF5 file.  Haplotypes should be stored in the slowest changing
#' dimension as defined in the HDF5 specification (note that different languages
#' treat this as rows or columns).  If the haplotypes are stored in the other
#' dimension then simply set the argument `transpose = TRUE`.
#' If the user is unsure of the convention of
#' the language they used to create the HDF5 file, then the simplest approach is
#' to simply load the data specifying only the HDF5 file name and then confirm
#' that number of haplotypes and their length have not been exchanged in the diagnostic
#' output which kalis prints.
#'
#' @param haps can be the name of a file from which the haplotypes are to be read, or
#'   can be a binary R matrix.
#' @param format the file format which `file` is stored in, or `"auto"` to
#'   detect format based on the file extension.  Recognised options are `"hdf5"`
#'   and `"vcf"`
#' @param ... options passed on to the specific format transcoding engine.  The following
#'   additional arguments are supported:
#'
#'   1. HDF5
#'       - `transpose` a logical value indicating whether to switch the
#'   interpretation of the slowest changing dimension (hence switching the
#'   number of haplotypes and the length)
#'
#' @return Nothing is returned by the function.  However, a status message is
#'   output indicating the dimensions of the loaded data.  This can be useful for
#'   HDF5 input when you are uncertain about the orientation of the `haps` object.
#'   If this message shows the dimensions incorrectly ordered,
#'   call the function again with the extra argument `transpose = TRUE`.
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
#' # When correct orientation is known, can avoid diagnostic messages for running
#' # in script files
#' suppressMessages(CacheHaplotypes("myhaps.h5"))
#' }
#'
#' @export
CacheHaplotypes <- function(haps, format = "auto", ...) {
  if(is.character(haps)) {
    if(!file.exists(haps)) {
      if(!is.na(get("N", envir = pkgVars))) {
        stop("Cannot find file to make new cache. Keeping existing cache.")
      } else {
        stop("Cannot find file.")
      }
    }

    cached <- FALSE
    ext <- stringr::str_to_lower(stringr::str_extract(haps, stringr::regex("\\.[0-9a-z]+$")))

    if(format == "hdf5" ||
       (format == "auto" && (ext == ".h5" ||
                             ext == ".hdf5"))) {
      CacheHaplotypes.hdf5(haps, ...)
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
    CacheHaplotypes.matrix(haps, ...)
  } else {
    stop("Unable to load haplotypes from the provided object.")
  }
}

CacheHaplotypes.matrix <- function(x, transpose = FALSE) {
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
  assign("N", N, envir = pkgVars)

  # Cache it!
  assign("L", CacheHaplotypes_matrix_2(x, N, L, transpose), envir = pkgVars)
}

CacheHaplotypes.hdf5 <- function(hdf5.file, transpose = FALSE) {
  # Check for file and dataset within file
  if(!file.exists(hdf5.file)) {
    if(!is.na(get("N", envir = pkgVars))) {
      stop("Cannot find HDF5 file to make new cache. Keeping existing cache.")
    } else {
      stop("Cannot find HDF5 file.")
    }
  }
  h5content <- h5ls(hdf5.file)
  if(!("haps" %in% h5content$name)) {
    if(!is.na(get("N", envir = pkgVars))) {
      stop("HDF5 file already exists but does not contain a 'haps' object in the root for haplotype data. Keeping existing cache.")
    } else {
      stop("HDF5 file already exists but does not contain a 'haps' object in the root for haplotype data.")
    }
  }

  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("N", envir = pkgVars))) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }


  # Get dimensions of haplotype data
  hdf5.dim <- as.integer(str_split_fixed(h5content[h5content$name=="haps","dim"], "x", n = Inf))
  # Sometimes on Linux it leaves a dimension unspecified, so have to load 1 row/col to ascertain
  if(is.na(hdf5.dim[1])) {
    hdf5.dim[1] <- dim(h5read(hdf5.file, "/haps", index = list(NULL,1)))[1]
  }
  if(is.na(hdf5.dim[2])) {
    hdf5.dim[2] <- dim(h5read(hdf5.file, "/haps", index = list(1,NULL)))[2]
  }
  if(!transpose) {
    N <- hdf5.dim[2]
    L <- hdf5.dim[1]
  } else {
    N <- hdf5.dim[1]
    L <- hdf5.dim[2]
  }

  # Hap IDs are only numeric for HDF5
  assign("N", N, envir = pkgVars)

  # Create a closure for easy access of HDF5 file
  # We'll read in 10MB (raw size) chunks
  make.hdf5.access <- function(hdf5.file, N, L) {
    step.size <- round(10000000/(L/8))
    current.step <- 1
    function() {
      if(current.step > N) {
        return(matrix(nrow = 0, ncol = 0))
      }
      upto <- min(current.step + step.size - 1, N)
      if(!transpose)
        res <- h5read(hdf5.file, "haps", index = list(NULL, current.step:upto))
      else
        res <- t(h5read(hdf5.file, "haps", index = list(current.step:upto, NULL)))
      current.step <<- upto + 1
      res
    }
  }

  # Cache it!
  assign("L", CacheHaplotypes_hdf5_2(make.hdf5.access(hdf5.file, N, L), N, L), envir = pkgVars)
}

#' Retrieve haplotypes from memory cache
#'
#' Retrieve haplotypes from the memory cache, converting the raw binary into
#' a simple integer 0/1 vector for use in R.
#'
#' To achieve higher performance, kalis internally represents haplotypes
#' in an efficient raw binary format in memory which cannot be directly viewed
#' or manipulated in R.  This function enables users to copy whole or partial
#' haplotypes out of this low-level format and into a standard R
#' vector of 0's and 1's.
#'
#' @param ids which haplotypes to retrieve from the cache
#' @param start the first locus position to retrieve for the specified haplotypes.
#'   Defaults to the beginning of the haplotype.
#' @param length the number of consecutive loci to return from the start position.
#'   Defaults to length of entire haplotype from `start` argument to end.
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
#' QueryCache(...)
#' }
#'
#' @export QueryCache
QueryCache <- function(ids = NA, start = 1, length = NA) {
  N <- get("N", envir = pkgVars)

  if(is.na(length)) {
    length <- get("L", envir = pkgVars)
  }
  if((length(ids)==1 && is.na(ids)) || all(ids==1:N)) {
    ids <- 1:N

    res <- matrix(nrow = length, ncol = length(ids))
    for(l in start:(start+length-1)) {
      res[l-start+1,] <- as.integer(intToBits(QueryCache2_loc(l-1)))[1:length(ids)]
    }
  } else {
    res <- matrix(nrow = length, ncol = length(ids))
    for(i in 1:length(ids)) {
      id <- ids[i]
      if(is.numeric(id) && abs(id - round(id)) < .Machine$double.eps^0.5) {
        if(id<1 || id>N) {
          stop("haplotype id ", id, " is out of range.")
        }
        idx <- id
      } else {
        stop(glue("ids must be integer values identifying haplotype number between 1 and {N}."))
      }

      hap <- as.integer(intToBits(QueryCache2_ind(idx-1)))

      res[,i] <- hap[start:(start+length-1)]
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
  ClearHaplotypeCache()
  library.dynam.unload("kalis", libpath)
}
