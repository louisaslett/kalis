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
