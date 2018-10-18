#' Write haplotype matrix to HDF5 formatted cache-friendly file
#'
#' Writes an R matrix of 0/1s to the HDF5 format which is used for reading to
#' optimised in memory cache.
#'
#' The primary method to load data into kalis' internal optimised cache is from
#' an HDF5 file storage.  If the user has a collection of haplotypes already
#' represented as a matrix of 0's and 1's in R, this function can be used to
#' write to HDF5 the format required to load into cache.
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
#' @export WriteIndividualSequenceH5
WriteIndividualSequenceH5 <- function(hdf5.file, ind.sequence) {
  if(any(!(ind.sequence == 1 | ind.sequence == 0))) {
    stop("sequence is not binary.")
  }
  if(!is.matrix(ind.sequence) && !is.vector(ind.sequence, mode = "numeric")) {
    stop("ind.sequence must be a vector (for one sequence) or a matrix (seq length x num seq).")
  }
  if(is.matrix(ind.sequence)) {
    N <- ncol(ind.sequence)
    L <- nrow(ind.sequence)
  } else {
    N <- 1
    L <- length(ind.sequence)
  }

  if(file.exists(hdf5.file)) {
    if(nrow(h5ls(hdf5.file) %>% filter(name == "seqs")) != 1) {
      stop("HDF5 file already exists but does not contain sequence data.")
    }
    hdf5.dim <- h5ls(hdf5.file) %>% filter(name == "seqs") %>% select(dim) %>% str_split_fixed("x", n = Inf) %>% as.integer()
    if(hdf5.dim[1] != L) {
      stop(glue("HDF5 file contains sequences of length {hdf5.dim[1]}, but ind.sequence contains sequences of length {L}."))
    }
    message("HDF5 file already exists, appending sequences ...\n")
  } else {
    message("Creating HDF5 file ...\n")
    h5createFile(hdf5.file)
    h5createDataset(file = hdf5.file,
                    dataset = "seqs",
                    dims = c(L, 0),
                    maxdims = c(L, 10e9),
                    storage.mode = "integer",
                    chunk = c(L, 1),
                    level = 7)
    hdf5.dim <- c(L, 0)
  }

  from <- hdf5.dim[2] + 1
  to <- hdf5.dim[2] + N

  # Expand to hold new sequences
  h5set_extent(hdf5.file, "seqs", c(hdf5.dim[1], to))

  # Write
  message(glue("Writing {N} sequence(s) of size {L} ...\n"))
  h5write(ind.sequence, hdf5.file, "seqs", index = list(NULL, from:to))
}

ReadIndividualSequenceH5 <- function(hdf5.file, inds) {
  if(!is.vector(inds, mode = "numeric")) {
    stop("inds must be a vector of sequence indices.")
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
  N <- hdf5.dim[2]
  L <- hdf5.dim[1]

  # Check index set range
  if(any(inds < 1) || any(inds > N)) {
    stop(glue("HDF5 file contains {N} sequences, some requested indices out of range."))
  }

  # Read
  ind.sequence <- h5read(hdf5.file, "seqs", index = list(NULL, inds))

  H5close()

  ind.sequence
}
