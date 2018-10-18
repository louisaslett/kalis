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
#' @export WriteIndividualHaplotypeH5
WriteIndividualHaplotypeH5 <- function(hdf5.file, ind.haplotype) {
  if(any(!(ind.haplotype == 1 | ind.haplotype == 0))) {
    stop("haplotype is not binary.")
  }
  if(!is.matrix(ind.haplotype) && !is.vector(ind.haplotype, mode = "numeric")) {
    stop("ind.haplotype must be a vector (for one haplotype) or a matrix (hap length x num hap).")
  }
  if(is.matrix(ind.haplotype)) {
    N <- ncol(ind.haplotype)
    L <- nrow(ind.haplotype)
  } else {
    N <- 1
    L <- length(ind.haplotype)
  }

  if(file.exists(hdf5.file)) {
    if(nrow(h5ls(hdf5.file) %>% filter(name == "haps")) != 1) {
      stop("HDF5 file already exists but does not contain haplotype data.")
    }
    hdf5.dim <- h5ls(hdf5.file) %>% filter(name == "haps") %>% select(dim) %>% str_split_fixed("x", n = Inf) %>% as.integer()
    if(hdf5.dim[1] != L) {
      stop(glue("HDF5 file contains haplotypes of length {hdf5.dim[1]}, but ind.haplotype contains haplotypes of length {L}."))
    }
    message("HDF5 file already exists, appending haplotypes ...\n")
  } else {
    message("Creating HDF5 file ...\n")
    h5createFile(hdf5.file)
    h5createDataset(file = hdf5.file,
                    dataset = "haps",
                    dims = c(L, 0),
                    maxdims = c(L, 10e9),
                    storage.mode = "integer",
                    chunk = c(L, 1),
                    level = 7)
    hdf5.dim <- c(L, 0)
  }

  from <- hdf5.dim[2] + 1
  to <- hdf5.dim[2] + N

  # Expand to hold new haplotypes
  h5set_extent(hdf5.file, "haps", c(hdf5.dim[1], to))

  # Write
  message(glue("Writing {N} haplotype(s) of size {L} ...\n"))
  h5write(ind.haplotype, hdf5.file, "haps", index = list(NULL, from:to))
}

ReadIndividualHaplotypeH5 <- function(hdf5.file, inds) {
  if(!is.vector(inds, mode = "numeric")) {
    stop("inds must be a vector of haplotype indices.")
  }

  # Check for file
  if(!file.exists(hdf5.file)) {
    stop("Cannot find HDF5 file.")
  }

  # Get dimensions of haplotype data
  if(nrow(h5ls(hdf5.file) %>% filter(name == "haps")) != 1) {
    stop("HDF5 file already exists but does not contain haplotype data.")
  }
  hdf5.dim <- h5ls(hdf5.file) %>% filter(name == "haps") %>% select(dim) %>% str_split_fixed("x", n = Inf) %>% as.integer()
  N <- hdf5.dim[2]
  L <- hdf5.dim[1]

  # Check index set range
  if(any(inds < 1) || any(inds > N)) {
    stop(glue("HDF5 file contains {N} haplotypes, some requested indices out of range."))
  }

  # Read
  ind.haplotype <- h5read(hdf5.file, "haps", index = list(NULL, inds))

  H5close()

  ind.haplotype
}
