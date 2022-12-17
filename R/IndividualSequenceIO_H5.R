#' I/O for haplotype matrices in HDF5 files
#'
#' Reads/writes an R matrix of 0/1s to the HDF5 format which is used for reading to the kalis optimised memory cache.
#' If you're working with a large haplotype dataset, we recommend that you convert it directly to this HDF5 format (see vignette) rather than read it into R.
#'
#' The primary method to load data into kalis' internal optimised cache is from an HDF5 storage file.
#' If the user has a collection of haplotypes already represented as a matrix of 0's and 1's in R, this function can be used to write to HDF5 in the format required to load into cache.
#'
#' kalis expects a 2-dimensional object named `haps` at the root level of the HDF5 file.
#' Haplotypes should be stored in the slowest changing dimension as defined in the HDF5 specification (note that different languages treat this as rows or columns).
#'
#' Note that if `hdf5.file` exists but does not contain a dataset named `haps`, then `WriteHaplotypes` will simply create a `haps` dataset within the existing file.
#'
#' @param hdf5.file the name of the file which the haplotypes are to be written to.
#' @param haps a vector or a matrix where each column is a haplotype to be stored in the file `hdf5.file`.
#' @param hap.ids a character vector naming haplotypes when writing, or which haplotypes are to be read.
#' @param loci.ids a character vector naming variants when writing, or which variants are to be read.
#' @param hap.idx an integer vector of the indices of which haplotypes are to be read (for naming, use `hap.ids`).
#' @param loci.idx an integer vector of the indices of which variants are to be read (for naming, use `hap.ids`).
#' @param haps.name a string providing the full path and object name where the haplotype matrix should be read/written.
#' @param hap.ids.name a string providing the full path and object name where the haplotype names (in `haps.ids`) should be read/written.
#' @param loci.ids.name a string providing the full path and object name where the variant names (in `loci.ids`) should be read/written.
#' @param append a logical indicating whether overwrite (default) or append to an existing `haps` dataset if it already exists in `hdf5.file`.
#' @param transpose a logical indicating whether to transpose the logic of haplotypes/variants when reading.
#'
#' @return `WriteHaplotypes` does not return anything.
#'
#'   `ReadHaplotypes` returns a binary matrix containing the
#' haplotypes that were specified in `ids`.
#'
#' @seealso [CacheHaplotypes()] to fill the kalis cache with haplotypes.
#'
#' @examples
#' \donttest{
#' # Generate a random mini set of haplotypes to write
#' n.haps <- 20
#' n.vars <- 200
#' haps <- matrix(sample(0:1, n.haps*n.vars, replace = TRUE),
#'                nrow = n.vars, ncol = n.haps)
#'
#' # ... write them to a file, giving alphabetic letters "A" through "T" as the #' # haplotype names ...
#' WriteHaplotypes("~/myhaps.h5", haps, hap.ids = LETTERS[1:20])
#'
#' # ... and confirm we can read a chosen portion back.  Try to read back
#' # the 10th and 11th haplotypes by using their name (J and K are 10th and 11th
#' # letter of the alphabet)
#' h5 <- ReadHaplotypes("~/myhaps.h5", hap.ids = c("J","K"))
#' all(h5$haps == haps[, 10:11])
#'
#' # Read from the .h5 file into the kalis cache and confirm that what we wrote
#' # out to the HDF5 file matches the original matrix we generated in R
#' CacheHaplotypes("~/myhaps.h5")
#' all(haps == QueryCache())
#' }
#'
#' @export WriteHaplotypes
WriteHaplotypes <- function(hdf5.file, haps,
                            hap.ids = NA, loci.ids = NA,
                            haps.name = "/haps", hap.ids.name = "/hap.ids", loci.ids.name = "/loci.ids",
                            append = FALSE) {
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
  N.old <- 0

  write.hap.ids <- FALSE
  if(!any(is.na(hap.ids))) {
    if(!is.atomic(hap.ids) || length(hap.ids) != N) {
      stop(glue("hap.ids must be a vector of length {N} to match supplied haplotype data."))
    }
    if(length(unique(hap.ids)) != length(hap.ids)) {
      stop("hap.ids must all be unique.")
    }
    write.hap.ids <- TRUE
  }

  write.loci.ids <- FALSE # see below
  if(!any(is.na(loci.ids))) {
    if(!is.atomic(loci.ids) || length(loci.ids) != L) {
      stop(glue("loci.ids must be a vector of length {L} to match supplied haplotype data."))
    }
    if(length(unique(loci.ids)) != length(loci.ids)) {
      stop("loci.ids must all be unique.")
    }
    # NB loci ids are special in that they are only written when writing a new
    #    dataset to disk or when overwriting a previous data set.  But when
    #    appending, we sanity check for matching but do not append (only append
    #    haps, not new loci!)
    write.loci.ids <- TRUE
  }

  if(file.exists(hdf5.file)) {
    h5 <- rhdf5::H5Fopen(hdf5.file)

    if(!rhdf5::H5Lexists(h5, haps.name)) {
      if(!any(is.na(hap.ids)) && rhdf5::H5Lexists(h5, hap.ids.name))
        stop(glue("HDF5 file exists and does not contain haplotype data, but already has a data set at {hap.ids.name} so cannot create haplotype IDs."))
      if(!any(is.na(loci.ids)) && rhdf5::H5Lexists(h5, loci.ids.name))
        stop(glue("HDF5 file exists and does not contain haplotype data, but already has a data set at {loci.ids.name} so cannot create locus IDs."))
      message(glue("HDF5 file already exists but does not contain haplotype data at {haps.name}, now adding."))
      if(!rhdf5::h5createDataset(file = hdf5.file,
                                 dataset = haps.name,
                                 dims = c(L, N),
                                 maxdims = c(L, rhdf5::H5Sunlimited()),
                                 H5type = "H5T_STD_U8LE",
                                 chunk = c(L, 1),
                                 level = 7))
        stop("Failed creating dataset.")
      if(!any(is.na(hap.ids))) {
        rhdf5::h5createDataset(file = hdf5.file,
                               dataset = hap.ids.name,
                               dims = c(N),
                               maxdims = c(rhdf5::H5Sunlimited()),
                               storage.mode = "character",
                               size = 100,
                               chunk = c(1000),
                               level = 7)
        write.hap.ids <- TRUE
      }
      if(!any(is.na(loci.ids))) {
        rhdf5::h5createDataset(file = hdf5.file,
                               dataset = loci.ids.name,
                               dims = c(L),
                               maxdims = c(rhdf5::H5Sunlimited()),
                               storage.mode = "character",
                               size = 100,
                               chunk = c(1000),
                               level = 7)
        write.loci.ids <- TRUE
      }
    } else {
      if(append) {
        message("HDF5 file already exists, appending haplotypes ...\n")

        if(!any(is.na(hap.ids)) && !rhdf5::H5Lexists(h5, hap.ids.name))
          stop(glue("There are no haplotype IDs in the {hap.ids.name} data set so cannot append haplotype IDs."))

        h5.haps <- rhdf5::H5Dopen(h5, haps.name)
        h5.haps.S <- rhdf5::H5Dget_space(h5.haps)
        haps.dims <- rhdf5::H5Sget_simple_extent_dims(h5.haps.S)$size
        if(length(haps.dims) != 2)
          stop(glue("The data set {haps.name} does not contain 2-dimensional data ... cannot append."))
        if(haps.dims[1] != L)
          stop(glue("The length of haplotypes in the {haps.name} data set is not {L} ... cannot append."))

        if(!any(is.na(hap.ids))) {
          if(!rhdf5::H5Lexists(h5, hap.ids.name)) {
            stop("haplotype IDs have been provided, but the existing haplotypes in the HDF5 file do not already have IDs.")
          }
          h5.hap.ids <- rhdf5::H5Dopen(h5, hap.ids.name)
          h5.hap.ids.S <- rhdf5::H5Dget_space(h5.hap.ids)
          hap.ids.dims <- rhdf5::H5Sget_simple_extent_dims(h5.hap.ids.S)$size
          if(length(hap.ids.dims) != 1)
            stop(glue("The data set {hap.ids.name} does not contain 1-dimensional data ... cannot append IDs."))
          if(hap.ids.dims[1] != N)
            stop(glue("The length of haplotype IDs in the {hap.ids.name} data set is not {N} ... cannot append IDs."))
          if(any(as.character(hap.ids) %in% as.character(rhdf5::h5read(h5, hap.ids.name)))) {
            stop("Some of the supplied hap.ids clash with IDs already in the data set.")
          }
          write.hap.ids <- TRUE
        }

        if(!any(is.na(loci.ids))) {
          if(!rhdf5::H5Lexists(h5, loci.ids.name)) {
            stop("locus IDs have been provided, but the existing loci in the HDF5 file do not already have IDs.")
          }
          h5.loci.ids <- rhdf5::H5Dopen(h5, loci.ids.name)
          h5.loci.ids.S <- rhdf5::H5Dget_space(h5.loci.ids)
          loci.ids.dims <- rhdf5::H5Sget_simple_extent_dims(h5.loci.ids.S)$size
          if(length(loci.ids.dims) != 1)
            stop(glue("The data set {loci.ids.name} does not contain 1-dimensional data ... invalid IDs."))
          if(loci.ids.dims[1] != L)
            stop(glue("The length of locus IDs in the {loci.ids.name} data set is not {L} ... invalid IDs."))
          h5.loci <- rhdf5::h5read(h5, loci.ids.name)
          if(any(as.character(loci.ids) != h5.loci))
            stop("The loci IDs provided and those in the HDF5 file differ.")
          write.loci.ids <- FALSE # we're appending *haps* so no new loci
        }

        # TODO: check data type for all these and error out if mismatch

        N.old <- haps.dims[2]

        if(write.hap.ids) {
          rhdf5::h5set_extent(h5, hap.ids.name, c(N.old+N))
        }
        rhdf5::h5set_extent(h5, haps.name, c(L, N.old+N))
      } else {
        message("HDF5 file exists and already contains a haps dataset, overwriting existing haps dataset...\n")

        # TODO: check data type for all these and error out if mismatch

        rhdf5::h5set_extent(h5, haps.name, c(L, N))
        if(!any(is.na(hap.ids)) && rhdf5::H5Lexists(h5, hap.ids.name))
          rhdf5::h5set_extent(h5, hap.ids.name, c(N))
        if(!any(is.na(loci.ids)) && rhdf5::H5Lexists(h5, loci.ids.name))
          rhdf5::h5set_extent(h5, loci.ids.name, c(L))
      }
    }
  } else {
    message("Creating HDF5 file ...\n")
    rhdf5::h5createFile(hdf5.file)

    if(!rhdf5::h5createDataset(file = hdf5.file,
                               dataset = haps.name,
                               dims = c(L, N),
                               maxdims = c(L, rhdf5::H5Sunlimited()),
                               H5type = "H5T_STD_U8LE",
                               chunk = c(L, 1),
                               level = 7))
      stop("Failed creating dataset.")
    if(!any(is.na(hap.ids)))
      rhdf5::h5createDataset(file = hdf5.file,
                             dataset = hap.ids.name,
                             dims = c(N),
                             maxdims = c(rhdf5::H5Sunlimited()),
                             storage.mode = "character",
                             size = 100,
                             chunk = c(1000),
                             level = 7)
    if(!any(is.na(loci.ids)))
      rhdf5::h5createDataset(file = hdf5.file,
                             dataset = loci.ids.name,
                             dims = c(L),
                             maxdims = c(rhdf5::H5Sunlimited()),
                             storage.mode = "character",
                             size = 100,
                             chunk = c(1000),
                             level = 7)

    h5 <- rhdf5::H5Fopen(hdf5.file)
  }

  # Write
  message(glue("Writing {N} haplotype(s) of size {L} ...\n"))

  rhdf5::h5write(as.array(haps), h5, haps.name, index = list(NULL, (N.old+1):(N.old+N)))
  if(write.hap.ids) {
    rhdf5::h5write(as.array(as.character(hap.ids)), h5, hap.ids.name, index = list((N.old+1):(N.old+N)))
  }
  if(write.loci.ids) {
    rhdf5::h5write(as.array(as.character(loci.ids)), h5, loci.ids.name, index = list(1:L))
  }

  rhdf5::h5closeAll()
}



IDs.to.index <- function(hap.ids, loci.ids, hap.idx, loci.idx, all.hap.ids, all.loci.ids) {
  if(!is.atomic(hap.ids) || !is.atomic(loci.ids) || !is.atomic(hap.idx) || !is.atomic(loci.idx)) {
    stop("Arguments can only be vectors.")
  }
  if(identical(hap.ids, NA) && identical(hap.idx, NA)) {
    stop("At least one of hap.ids or hap.idx must be provided.")
  }
  if(!identical(hap.ids, NA) && !identical(hap.idx, NA)) {
    stop("Only one of hap.ids or hap.idx may be provided.")
  }
  if(identical(loci.ids, NA) && identical(loci.idx, NA)) {
    stop("At least one of loci.ids or loci.idx must be provided.")
  }
  if(!identical(loci.ids, NA) && !identical(loci.idx, NA)) {
    stop("Only one of loci.ids or loci.idx may be provided.")
  }

  # Check IDs and argument compatibility
  if(!identical(hap.ids, NA) && is.integer(all.hap.ids)) {
    hap.ids <- suppressWarnings(as.integer(hap.ids))
    if(any(is.na(hap.ids))) {
      stop("hap.ids supplied but no hap.ids are present ... failed when trying to interpret as an index.")
    }
    warning("hap.ids supplied but no hap.ids are present ... interpreting as an index.")
  } else if(!identical(hap.ids, NA) && is.character(all.hap.ids)) {
    hap.ids2 <- match(as.character(hap.ids), all.hap.ids)
    if(any(is.na(hap.ids2))) {
      stop(glue("Can't find haplotype IDs: {paste(as.character(hap.ids)[is.na(hap.ids2)], collapse = ', ')}"))
    }
    hap.ids <- hap.ids2
  } else if(identical(hap.ids, NA)) {
    hap.ids <- suppressWarnings(as.integer(hap.idx))
    if(any(is.na(hap.ids))) {
      stop("Failed when trying to interpret hap.idx as an integer.")
    }
  } else {
    stop("Unrecoverable error trying to interpret hap.ids/hap.idx arguments.")
  }

  if(!identical(loci.ids, NA) && is.integer(all.loci.ids)) {
    loci.ids <- suppressWarnings(as.integer(loci.ids))
    if(any(is.na(loci.ids))) {
      stop("loci.ids supplied but no loci.ids are present ... failed when trying to interpret as an index.")
    }
    warning("loci.ids supplied but no loci.ids are present ... interpreting as an index.")
  } else if(!identical(loci.ids, NA) && is.character(all.loci.ids)) {
    loci.ids2 <- match(as.character(loci.ids), all.loci.ids)
    if(any(is.na(loci.ids2))) {
      stop(glue("loci IDs: {paste(as.character(loci.ids)[is.na(loci.ids2)], collapse = ', ')} not found."))
    }
    loci.ids <- loci.ids2
  } else if(identical(loci.ids, NA)) {
    loci.ids <- suppressWarnings(as.integer(loci.idx))
    if(any(is.na(loci.ids))) {
      stop("Failed when trying to interpret loci.idx as an integer.")
    }
  } else {
    stop("Unrecoverable error trying to interpret hap.ids/hap.idx arguments.")
  }

  # Check we have sensible indexing by this point
  if(any(hap.ids < 1 | hap.ids > length(all.hap.ids))) {
    stop("Invalid haplotypes specified.")
  }
  if(any(loci.ids < 1 | loci.ids > length(all.loci.ids))) {
    stop("Invalid loci specified.")
  }

  list(hap = hap.ids, loci = loci.ids)
}



#' @rdname WriteHaplotypes
#' @export ReadHaplotypes
ReadHaplotypes <- function(hdf5.file,
                           loci.idx = NA, hap.idx = NA,
                           loci.ids = NA, hap.ids = NA,
                           haps.name = "/haps",
                           loci.ids.name = "/loci.ids", hap.ids.name = "/hap.ids",
                           transpose = FALSE) {
  if(!file.exists(hdf5.file)) {
    stop("Cannot find HDF5 file.")
  }
  if(!identical(hap.ids, NA) && !identical(hap.idx, NA)) {
    stop("Can only specify one of hap.ids or hap.idx argument.")
  }
  if(!identical(loci.ids, NA) && !identical(loci.idx, NA)) {
    stop("Can only specify one of loci.ids or loci.idx argument.")
  }

  # Get dimensions of hap data
  hdf5.dim <- integer(2)
  hdf5.dim[1] <- dim(rhdf5::h5read(hdf5.file, haps.name, index = list(NULL,1)))[1]
  hdf5.dim[2] <- dim(rhdf5::h5read(hdf5.file, haps.name, index = list(1,NULL)))[2]
  if(!transpose) {
    N <- hdf5.dim[2]
    L <- hdf5.dim[1]
  } else {
    N <- hdf5.dim[1]
    L <- hdf5.dim[2]
  }

  # Figure out what to call all haplotypes/loci
  all.hap.ids  <- tryCatch(rhdf5::h5read(hdf5.file, hap.ids.name),  error = function(e) 1:N)
  all.loci.ids <- tryCatch(rhdf5::h5read(hdf5.file, loci.ids.name), error = function(e) 1:L)

  # Default to get everything
  if(identical(hap.ids, NA) && identical(hap.idx, NA)) {
    hap.idx <- 1:N
  }
  if(identical(loci.ids, NA) && identical(loci.idx, NA)) {
    loci.idx <- 1:L
  }
  idx <- IDs.to.index(hap.ids, loci.ids, hap.idx, loci.idx, all.hap.ids, all.loci.ids)

  # Read
  if(!transpose)
    haps <- matrix(as.integer(rhdf5::h5read(hdf5.file, haps.name, index = list(idx$loci, idx$hap))), nrow = length(idx$loci))
  else
    haps <- matrix(as.integer(t(rhdf5::h5read(hdf5.file, haps.name, index = list(idx$hap, idx$loci)))), nrow = length(idx$loci))

  rhdf5::H5close()

  list(haps = haps, hap.ids = all.hap.ids[idx$hap], loci.ids = all.loci.ids[idx$loci])
}
