#' Write haplotype matrix to HDF5 file
#'
#' Writes an R matrix of 0/1s to the HDF5 format which is used for reading to
#' optimised in memory cache.  If you're working with a large haplotype dataset,
#' we recommend that you convert it directly to this HDF5 format (see vignette)
#' rather than read it into R.
#'
#' The primary method to load data into kalis' internal optimised cache is from
#' an HDF5 storage file.  If the user has a collection of haplotypes already
#' represented as a matrix of 0's and 1's in R, this function can be used to
#' write to HDF5 in the format required to load into cache.
#'
#' kalis expects a 2-dimensional object named \code{haps} at the root
#' level of the HDF5 file.  Haplotypes should be stored in the slowest changing
#' dimension as defined in the HDF5 specification (note that different languages
#' treat this as rows or columns).
#'
#' Note that if \code{hdf5.file} exists but does not contain a dataset named
#' \code{haps}, then \code{WriteIndividualHaplotypeH5} will simply create a
#' \code{haps} dataset within the existing file.
#'
#' @param hdf5.file the name of the file which the haplotypes are to be
#'   written to.
#' @param haps a vector or a matrix where each column is a haplotype to be
#'   stored in the file \code{hdf5.file}.
#' @param ids a vector of the indices of which haplotypes are to be read.
#' @param append a logical indicating whether overwrite (default) or append to
#'   an existing \code{haps} dataset if it already exists in \code{hdf5.file}.
#'
#' @return
#' \code{WriteIndividualHaplotypeH5} does not return anything.
#'
#' \code{ReadIndividualHaplotypeH5} returns a binary matrix containing the
#' haplotypes that were specified in \code{ids}.
#'
#' @seealso
#' \code{\link{CacheHaplotypes}} to fill the kalis cache with haplotypes
#'
#' @examples
#' # Examples
#' \dontrun{
#' # For the purposes of an example, generate random haplotypes ...
#' n.haps <- 100
#' n.loci <- 20000
#' haps <- matrix(sample(0:1, n.haps*n.loci, replace = TRUE),
#'                nrow = n.loci,
#'                ncol = n.haps)
#'
#' # ... write them to a file ...
#' WriteIndividualHaplotypeH5("myhaps.h5", haps)
#'
#' # ... and confirm we can read a chosen portion back.  Try to read back
#' # the 10th and 11th haplotypes
#' res <- ReadIndividualHaplotypeH5("myhaps.h5", 10:11)
#' all(res == haps[, 10:11])
#' }
#'
#' @export WriteIndividualHaplotypeH5
WriteIndividualHaplotypeH5 <- function(hdf5.file, haps, append = FALSE) {
  if(any(!(haps == 1 | haps == 0))) {
    stop("haplotype is not binary.")
  }
  if(!is.matrix(haps) && !is.vector(haps, mode = "numeric")) {
    stop("haps must be a vector (for one haplotype) or a matrix (hap length x num hap).")
  }
  if(is.matrix(haps)) {
    N <- ncol(haps)
    L <- nrow(haps)
  } else {
    N <- 1
    L <- length(haps)
  }

  if(file.exists(hdf5.file)) {
    h5content <- rhdf5::h5ls(hdf5.file)
    if(!("haps" %in% h5content$name)) {
      message("HDF5 file already exists but does not contain haplotype data, now adding haps dataset.")
      rhdf5::h5createDataset(file = hdf5.file,
                             dataset = "haps",
                             dims = c(L, 0),
                             maxdims = c(L, 10e9),
                             H5type = "H5T_STD_U8LE",
                             chunk = c(L, 1),
                             level = 7)
      hdf5.dim <- c(L, 0)
    } else {
      if(append) {
        hdf5.dim <- as.integer(str_split_fixed(h5content[h5content$name=="haps","dim"], "x", n = Inf))
        if(hdf5.dim[1] != L) {
          stop(glue("HDF5 file contains haplotypes of length {hdf5.dim[1]}, but haps contains haplotypes of length {L}."))
        }
        message("HDF5 file already exists, appending haplotypes ...\n")
      } else {
        message("HDF5 file exists and already contains a haps dataset, overwriting existing haps dataset...\n")
        rhdf5::h5delete(file = hdf5.file, name = "haps")
        rhdf5::h5createDataset(file = hdf5.file,
                               dataset = "haps",
                               dims = c(L, 0),
                               maxdims = c(L, 10e9),
                               H5type = "H5T_STD_U8LE",
                               chunk = c(L, 1),
                               level = 7)
        hdf5.dim <- c(L, 0)
      }

    }

  } else {
    message("Creating HDF5 file ...\n")
    rhdf5::h5createFile(hdf5.file)
    rhdf5::h5createDataset(file = hdf5.file,
                           dataset = "haps",
                           dims = c(L, 0),
                           maxdims = c(L, 10e9),
                           H5type = "H5T_STD_U8LE",
                           chunk = c(L, 1),
                           level = 7)
    hdf5.dim <- c(L, 0)
  }

  from <- hdf5.dim[2] + 1
  to <- hdf5.dim[2] + N

  # Expand to hold new haplotypes
  rhdf5::h5set_extent(hdf5.file, "haps", c(hdf5.dim[1], to))

  # Write
  message(glue("Writing {N} haplotype(s) of size {L} ...\n"))
  rhdf5::h5write(haps, hdf5.file, "haps", index = list(NULL, from:to))
}


#' @describeIn WriteIndividualHaplotypeH5 Read haplotype matrix from HDF5 file
#' @export ReadIndividualHaplotypeH5
ReadIndividualHaplotypeH5 <- function(hdf5.file, ids) {
  if(!is.vector(ids, mode = "numeric")) {
    stop("ids must be a vector of haplotype indices.")
  }

  # Check for file
  if(!file.exists(hdf5.file)) {
    stop("Cannot find HDF5 file.")
  }

  # Get dimensions of haplotype data
  h5content <- rhdf5::h5ls(hdf5.file)
  if(!("haps" %in% h5content$name)) {
    stop("HDF5 file already exists but does not contain a 'haps' object in the root for haplotype data.")
  }
  hdf5.dim <- as.integer(str_split_fixed(h5content[h5content$name=="haps","dim"], "x", n = Inf))
  N <- hdf5.dim[2]
  L <- hdf5.dim[1]

  # Check index set range
  if(any(ids < 1) || any(ids > N)) {
    stop(glue("HDF5 file contains {N} haplotypes, some requested indices out of range."))
  }

  # Read
  haps <- rhdf5::h5read(hdf5.file, "haps", index = list(NULL, ids))

  rhdf5::H5close()

  haps
}
