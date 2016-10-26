WriteIndividualSequence <- function(ind.sequence, id = NA) {
  working.dir <- get("working.dir", envir = pkgCache)

  if(any(!(ind.sequence == 1 | ind.sequence == 0))) {
    stop("sequence is not binary.")
  }

  if(is.na(id)) {
    id <- paste(sample(c(0:9, letters, LETTERS), 8, replace=TRUE), collapse = "")
    while(file.exists(paste0(working.dir, "/", id, ".seq"))) {
      id <- paste(sample(c(0:9, letters, LETTERS), 8, replace=TRUE), collapse = "")
    }
  }
  if(file.exists(paste0(working.dir, "/", id, ".seq"))) {
    warning("individual ", id, " already exists ... overwriting sequence data.")
  }

  #system(paste0("touch ", working.dir, "/", id, ".seq"))
  fd <- file(paste0(working.dir, "/", id, ".seq"), "wb")

  # Record the sequence length at start of the file
  writeBin(as.integer(length(ind.sequence)), fd, endian = "little")

  # Correct length of sequence to be byte sized
  if(length(ind.sequence) %% 8 != 0) {
    ind.sequence <- c(ind.sequence, rep(0, 8 - length(ind.sequence) %% 8))
  }

  # Write out to disk
  writeBin(packBits(as.logical(ind.sequence)), fd, endian = "little")

  close(fd)
}

ReadIndividualSequence <- function(id) {
  working.dir <- get("working.dir", envir = pkgCache)

  if(!file.exists(paste0(working.dir, "/", id, ".seq"))) {
    stop("individual ", id, " does not exist in the current working directory")
  }

  fd <- file(paste0(working.dir, "/", id, ".seq"), "rb")

  # Find out sequence length
  n <- readBin(fd, integer(), 1, endian = "little")

  # Extract sequence
  ind.sequence <- as.integer(rawToBits(readBin(fd, raw(), ceiling(n/8), endian = "little")))[1:n]

  close(fd)

  ind.sequence
}
