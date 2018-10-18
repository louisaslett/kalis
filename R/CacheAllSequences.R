pkgCache <- new.env()
assign("working.dir", '.', envir = pkgCache)
assign("seqs", NA, envir = pkgCache)
assign("seq_size", NA, envir = pkgCache)

CacheAllSequences <- function() {
  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("seqs", envir = pkgCache)[1])) {
    warning("sequences already cached ... overwriting existing cache.")
    ClearSequenceCache()
  }

  working.dir <- get("working.dir", envir = pkgCache)

  # Save index of sequence ids
  assign("seqs", sub(".seq", "", list.files(GetGenWorkingDir(), "*.seq")), envir = pkgCache)

  # Get all the file locations for the sequences and their size
  seqs <- list.files(GetGenWorkingDir(), "*.seq", full.names = TRUE)
  sizes <- file.size(seqs)
  if(!all(sizes == sizes[1])) {
    error("not all sequences in the working directory are the same length.")
  }
  buf.size <- sum(file.size(seqs)) # strictly speaking, could remove all the space at front for ints

  # Cache it!
  assign("seq_size", CacheAllSequences2(seqs, buf.size), envir = pkgCache)
}

# We assume haplotypes are stored in the slowest changing dimension
# per the HDF5 spec definition.  This is "row-wise" in the C standard
# spec, or "col-wise" in the rhdf5 spec.

#' Create memory cache of haplotype sequences
#'
#' Load haplotype sequences from hard drive to memory.
#'
#' To achieve higher performance, kalis internally represents haplotype sequences
#' in an efficient raw binary format in memory.  This function will load sequences
#' from an HDF5 file and convert this into kalis' internal format ready for use
#' by the other functions in this package.
#'
#' The package expects a 2-dimensional object named \code{seqs} at the root
#' level of the HDF5 file.  Haplotypes should be stored in the slowest changing
#' dimension as defined in the HDF5 specification (note that different languages
#' treat this as rows or columns).  If the haplotypes are stored in the other
#' dimension then simply set the argument \code{transpose = TRUE}.
#' If the user is unsure of the convention of
#' the language they used to create the HDF5 file, then the simplest approach is
#' to just load the data specifying only the HDF5 file name and then confirm
#' that number of haplotypes and sequence length have not been exchanged.
#'
#' @param from_recipient first recipient haplotype if creating a partial forward
#'   table.  By default includes from the first recipient haplotype.
#' @param to_recipient last recipient haplotype if creating a partial forward
#'   table.  By default includes to the last recipient haplotype.
#'
#' @return Nothing is returned by the function.  However, a status message is
#'   output indicating how the dimensions of the HDF5 \code{seqs} object have
#'   been interpreted.  If this message shows the dimensions incorrectly ordered,
#'   call the function again with argument \code{transpose = TRUE}.
#'
#' @seealso \code{\link{Forward}} to propagate the newly created table forward
#'   through the genome.
#'
#' @examples
#' # Examples
#' \dontrun{
#'
#' # Load sequences to cache from an HDF5 file on disk
#' CacheAllSequencesH5("myseqs.h5")
#'
#' # If the diagnostic message printed during the above indicates the numbers
#' # of haplotypes and sequence length are the wrong way around, reload with
#' # argument to transpose
#' CacheAllSequencesH5("myseqs.h5", transpose = TRUE)
#'
#' # When correct orientation is known, can avoid diagnostic messages for running
#' # in script files
#' suppressMessages(CacheAllSequencesH5("myseqs.h5"))
#' }
#'
#' @export CacheAllSequencesH5
CacheAllSequencesH5 <- function(hdf5.file, transpose = FALSE) {
  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("seqs", envir = pkgCache)[1])) {
    warning("sequences already cached ... overwriting existing cache.")
    ClearSequenceCache()
  }

  # Check for file
  if(!file.exists(hdf5.file)) {
    stop("Cannot find HDF5 file.")
  }

  # Get dimensions of sequence data
  if(nrow(h5ls(hdf5.file) %>% filter(name == "seqs")) != 1) {
    stop("HDF5 file already exists but does not contain sequence data.")
  }
  hdf5.dim <- h5ls(hdf5.file) %>% filter(name == "seqs") %>% select(dim) %>% str_split_fixed("x", n = Inf) %>% as.integer()
  if(!transpose) {
    N <- hdf5.dim[2]
    L <- hdf5.dim[1]
  } else {
    N <- hdf5.dim[1]
    L <- hdf5.dim[2]
  }

  # Seq IDs are only numeric for HDF5
  assign("seqs", 1:N, envir = pkgCache)

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
        res <- h5read(hdf5.file, "seqs", index = list(NULL, current.step:upto))
      else
        res <- t(h5read(hdf5.file, "seqs", index = list(current.step:upto, NULL)))
      current.step <<- upto + 1
      res
    }
  }

  # Cache it!
  assign("seq_size", CacheAllSequencesH52(make.hdf5.access(hdf5.file, N, L), N, L), envir = pkgCache)
}

#' Retrieve haplotypes from memory cache
#'
#' Retrieve haplotypes from the memory cache, converting the raw binary into
#' a simple integer 0/1 vector for use in R.
#'
#' To achieve higher performance, kalis internally represents haplotype sequences
#' in an efficient raw binary format in memory which cannot be directly viewed
#' or manipulated in R.  This function enables users to copy whole or partial
#' haplotype sequences out of this low-level format and into a standard R
#' vector of 0's and 1's.
#'
#' @param ids which haplotypes to retrieve from the cache
#' @param start the first locus position to retrieve for the specified haplotypes
#' @param length the number of consecutive loci to return from the start position
#'
#' @return A matrix of 0/1 integers with \code{length} rows and
#'   \code{length(ids)} columns, such that haplotypes appear in columns.
#'
#' @seealso \code{\link{Forward}} to propagate the newly created table forward
#'   through the genome.
#'
#' @examples
#' # Examples
#' \dontrun{
#' }
#'
#' @export QueryCache
QueryCache <- function(ids = NA, start = 1, length = NA) {
  seqs <- get("seqs", envir = pkgCache)

  if(is.na(length)) {
    length <- get("seq_size", envir = pkgCache)
  }
  if((length(ids)==1 && is.na(ids)) || all(ids==1:length(seqs))) {
    ids <- 1:length(seqs)

    res <- matrix(nrow = length(ids), ncol = length)
    for(l in start:(start+length-1)) {
      res[,l-start+1] <- as.integer(intToBits(QueryCache2_loc(l-1)))[1:length(ids)]
    }
  } else {
    res <- matrix(nrow = length(ids), ncol = length)
    for(i in 1:length(ids)) {
      id <- ids[i]
      if(is.numeric(id) && abs(id - round(id)) < .Machine$double.eps^0.5) {
        if(id<1 || id>length(seqs)) {
          stop("sequence id ", id, " is out of range.")
        }
        idx <- id
      } else {
        idx <- which(seqs == id)
        if(length(idx)==0) {
          stop("sequence id ", id, " is not in the cache.")
        }
      }

      seq <- as.integer(intToBits(QueryCache2_ind(idx-1)))

      res[i,] <- seq[start:(start+length-1)]
    }
  }

  res
}

#' Remove all cached sequences and free memory
#'
#' Remove all the sequences which were cached using a previous caching call and
#' free the memory that was allocated up for future use.
#'
#' To achieve higher performance, kalis internally represents haplotype sequences
#' in an efficient raw binary format in memory which cannot be directly viewed
#' or manipulated in R.  In particular, this cache sits outside R's memory
#' management and will never be garbage collected (unless R is quit or the
#' package is unloaded).  Therefore, this function is provided to enable freeing
#' the memory used by this cache.
#'
#' @return Nothing is returned.
#'
#' @seealso \code{\link{Forward}} to propagate the newly created table forward
#'   through the genome.
#'
#' @examples
#' # Examples
#' \dontrun{
#' }
#'
#' @export ClearSequenceCache
ClearSequenceCache <- function() {
  assign("seqs", NA, envir = pkgCache)
  assign("seq_size", NA, envir = pkgCache)
  ClearSequenceCache2()
}

.onUnload <- function(libpath) {
  ClearSequenceCache()
  library.dynam.unload("kalis", libpath)
}
