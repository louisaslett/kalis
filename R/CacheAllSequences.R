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

ClearSequenceCache <- function() {
  assign("seqs", NA, envir = pkgCache)
  assign("seq_size", NA, envir = pkgCache)
  ClearSequenceCache2()
}

.onUnload <- function(libpath) {
  ClearSequenceCache()
  library.dynam.unload("kalis", libpath)
}
