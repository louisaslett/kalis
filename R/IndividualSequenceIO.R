WriteIndividualHaplotype <- function(ind.haplotype, id = NA) {
  working.dir <- get("working.dir", envir = pkgCache)

  if(any(!(ind.haplotype == 1 | ind.haplotype == 0))) {
    stop("haplotype is not binary.")
  }

  if(is.na(id)) {
    id <- paste(sample(c(0:9, letters, LETTERS), 8, replace=TRUE), collapse = "")
    while(file.exists(paste0(working.dir, "/", id, ".hap"))) {
      id <- paste(sample(c(0:9, letters, LETTERS), 8, replace=TRUE), collapse = "")
    }
  }
  if(file.exists(paste0(working.dir, "/", id, ".hap"))) {
    warning("individual ", id, " already exists ... overwriting haplotype data.")
  }

  #system(paste0("touch ", working.dir, "/", id, ".hap"))
  fd <- file(paste0(working.dir, "/", id, ".hap"), "wb")

  # Record the haplotype length at start of the file
  writeBin(as.integer(length(ind.haplotype)), fd, endian = "little")

  # Correct length of haplotype to be byte sized
  if(length(ind.haplotype) %% 32 != 0) {
    ind.haplotype <- c(ind.haplotype, rep(0, 32 - length(ind.haplotype) %% 32))
  }

  # Write out to disk
  writeBin(packBits(as.logical(ind.haplotype), type = "integer"), fd, endian = "little")

  close(fd)
}

ReadIndividualHaplotype <- function(id) {
  working.dir <- get("working.dir", envir = pkgCache)

  if(!file.exists(paste0(working.dir, "/", id, ".hap"))) {
    stop("individual ", id, " does not exist in the current working directory")
  }

  fd <- file(paste0(working.dir, "/", id, ".hap"), "rb")

  # Find out haplotype length
  n <- readBin(fd, integer(), 1, endian = "little")

  # Extract haplotype
  ind.haplotype <- as.integer(intToBits(readBin(fd, "int", ceiling(n/32), endian = "little")))[1:n]

  close(fd)

  ind.haplotype
}
