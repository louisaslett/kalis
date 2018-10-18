pkgCache <- new.env()
assign("working.dir", '.', envir = pkgCache)
assign("haps", NA, envir = pkgCache)
assign("hap_size", NA, envir = pkgCache)

CacheAllHaplotypes <- function() {
  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("haps", envir = pkgCache)[1])) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }

  working.dir <- get("working.dir", envir = pkgCache)

  # Save index of haplotype ids
  assign("haps", sub(".hap", "", list.files(GetGenWorkingDir(), "*.hap")), envir = pkgCache)

  # Get all the file locations for the haplotypes and their size
  haps <- list.files(GetGenWorkingDir(), "*.hap", full.names = TRUE)
  sizes <- file.size(haps)
  if(!all(sizes == sizes[1])) {
    error("not all haplotypes in the working directory are the same length.")
  }
  buf.size <- sum(file.size(haps)) # strictly speaking, could remove all the space at front for ints

  # Cache it!
  assign("hap_size", CacheAllHaplotypes2(haps, buf.size), envir = pkgCache)
}

# We assume haplotypes are stored in the slowest changing dimension
# per the HDF5 spec definition.  This is "row-wise" in the C standard
# spec, or "col-wise" in the rhdf5 spec.

#' Create memory cache of haplotypes
#'
#' Load haplotypes from hard drive to memory.
#'
#' To achieve higher performance, kalis internally represents haplotypes
#' in an efficient raw binary format in memory.  This function will load haplotypes
#' from an HDF5 file and convert this into kalis' internal format ready for use
#' by the other functions in this package.
#'
#' The package expects a 2-dimensional object named \code{haps} at the root
#' level of the HDF5 file.  Haplotypes should be stored in the slowest changing
#' dimension as defined in the HDF5 specification (note that different languages
#' treat this as rows or columns).  If the haplotypes are stored in the other
#' dimension then simply set the argument \code{transpose = TRUE}.
#' If the user is unsure of the convention of
#' the language they used to create the HDF5 file, then the simplest approach is
#' to just load the data specifying only the HDF5 file name and then confirm
#' that number of haplotypes and their length have not been exchanged.
#'
#' @param hdf5.file the name of the file which the haplotypes are to be read from.
#' @param transpose a logical value indicating whether to switch the
#'   interpretation of the slowest changing dimension (hence switching the
#'   number of haplotypes and the length)
#'
#' @return Nothing is returned by the function.  However, a status message is
#'   output indicating how the dimensions of the HDF5 \code{haps} object have
#'   been interpreted.  If this message shows the dimensions incorrectly ordered,
#'   call the function again with argument \code{transpose = TRUE}.
#'
#' @seealso \code{\link{ClearHaplotypeCache}} to remove the haplotypes from
#'   the cache and free the memory.
#'
#' @examples
#' # Examples
#' \dontrun{
#'
#' # Load haplotypes to cache from an HDF5 file on disk
#' CacheAllHaplotypesH5("myhaps.h5")
#'
#' # If the diagnostic message printed during the above indicates the numbers
#' # of haplotypes and their length are the wrong way around, reload with
#' # argument to transpose
#' CacheAllHaplotypesH5("myhaps.h5", transpose = TRUE)
#'
#' # When correct orientation is known, can avoid diagnostic messages for running
#' # in script files
#' suppressMessages(CacheAllHaplotypesH5("myhaps.h5"))
#' }
#'
#' @export CacheAllHaplotypesH5
CacheAllHaplotypesH5 <- function(hdf5.file, transpose = FALSE) {
  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("haps", envir = pkgCache)[1])) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }

  # Check for file
  if(!file.exists(hdf5.file)) {
    stop("Cannot find HDF5 file.")
  }

  # Get dimensions of haplotype data
  if(nrow(h5ls(hdf5.file) %>% filter(name == "haps")) != 1) {
    stop("HDF5 file already exists but does not contain a 'haps' object in the root for haplotype data.")
  }
  hdf5.dim <- h5ls(hdf5.file) %>% filter(name == "haps") %>% select(dim) %>% str_split_fixed("x", n = Inf) %>% as.integer()
  if(!transpose) {
    N <- hdf5.dim[2]
    L <- hdf5.dim[1]
  } else {
    N <- hdf5.dim[1]
    L <- hdf5.dim[2]
  }

  # Hap IDs are only numeric for HDF5
  assign("haps", 1:N, envir = pkgCache)

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
  assign("hap_size", CacheAllHaplotypesH52(make.hdf5.access(hdf5.file, N, L), N, L), envir = pkgCache)
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
#'   Defaults to length of entire haplotype from \code{start} argument to end.
#'
#' @return A matrix of 0/1 integers with \code{length} rows and
#'   \code{length(ids)} columns, such that haplotypes appear in columns.
#'
#' @seealso \code{\link{CacheAllHaplotypesH5}} to fill the memory cache with
#'   haplotypes.
#'
#' @examples
#' # Examples
#' \dontrun{
#' }
#'
#' @export QueryCache
QueryCache <- function(ids = NA, start = 1, length = NA) {
  haps <- get("haps", envir = pkgCache)

  if(is.na(length)) {
    length <- get("hap_size", envir = pkgCache)
  }
  if((length(ids)==1 && is.na(ids)) || all(ids==1:length(haps))) {
    ids <- 1:length(haps)

    res <- matrix(nrow = length(ids), ncol = length)
    for(l in start:(start+length-1)) {
      res[,l-start+1] <- as.integer(intToBits(QueryCache2_loc(l-1)))[1:length(ids)]
    }
  } else {
    res <- matrix(nrow = length(ids), ncol = length)
    for(i in 1:length(ids)) {
      id <- ids[i]
      if(is.numeric(id) && abs(id - round(id)) < .Machine$double.eps^0.5) {
        if(id<1 || id>length(haps)) {
          stop("haplotype id ", id, " is out of range.")
        }
        idx <- id
      } else {
        idx <- which(haps == id)
        if(length(idx)==0) {
          stop("haplotype id ", id, " is not in the cache.")
        }
      }

      hap <- as.integer(intToBits(QueryCache2_ind(idx-1)))

      res[i,] <- hap[start:(start+length-1)]
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
#' @seealso \code{\link{CacheAllHaplotypesH5}} to propagate the newly created table forward
#'   through the genome.
#'
#' @examples
#' # Examples
#' \dontrun{
#' }
#'
#' @export ClearHaplotypeCache
ClearHaplotypeCache <- function() {
  assign("haps", NA, envir = pkgCache)
  assign("hap_size", NA, envir = pkgCache)
  ClearHaplotypeCache2()
}

.onUnload <- function(libpath) {
  ClearHaplotypeCache()
  library.dynam.unload("kalis", libpath)
}
