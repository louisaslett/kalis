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

QueryCache <- function(id, start = 1, length = NA) {
  seqs <- get("seqs", envir = pkgCache)

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

  seq <- as.integer(rawToBits(QueryCache2(idx-1)))
  if(is.na(length)) {
    length <- get("seq_size", envir = pkgCache)
  }

  seq[start:(start+length-1)]
}

ClearSequenceCache <- function() {
  assign("seqs", NA, envir = pkgCache)
  assign("seq_size", NA, envir = pkgCache)
  ClearSequenceCache2()
}

.onUnload <- function(libpath) {
  ClearSequenceCache()
  library.dynam.unload("StatGen", libpath)
}
