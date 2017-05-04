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
  library.dynam.unload("StatGen", libpath)
}
