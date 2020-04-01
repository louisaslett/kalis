CacheHaplotypes.hdf5 <- function(hdf5.file,
                                 loci.idx,
                                 hap.idx,
                                 transpose = FALSE,
                                 haps = "/haps",
                                 hdf5.pkg = c("hdf5r", "rhdf5")) {
  for(pkg in hdf5.pkg) {
    if(!(pkg %in% c("hdf5r", "rhdf5"))) {
      stop(glue("The package {pkg} is not supported for reading from HDF5 files."))
    }
    if(length(find.package(pkg, quiet = TRUE)) > 0) {
      if(pkg == "hdf5r") {
        return(CacheHaplotypes.hdf5.hdf5r(hdf5.file, loci.idx, hap.idx, transpose, haps))
      }
      if(pkg == "rhdf5") {
        return(CacheHaplotypes.hdf5.rhdf5(hdf5.file, loci.idx, hap.idx, transpose, haps))
      }
    }
  }
}



CacheHaplotypes.hdf5.hdf5r <- function(hdf5.file,
                                       loci.idx,
                                       hap.idx,
                                       transpose,
                                       haps) {
  # Check for file and dataset within file
  if(!testFile(hdf5.file, access = "r")) {
    CacheHaplotypes.err(checkFile(hdf5.file, access = "r"))
  }

  # Open file and objects
  h5 <- hdf5r::H5File$new(hdf5.file, "r")
  if(!h5$exists(haps)) {
    h5$close_all()
    CacheHaplotypes.err(glue("HDF5 file exists but does not contain an object at '{haps}' for the haplotype data."))
  }

  # locate haplotypes within file
  h5.haps <- h5$open(haps)
  if(length(h5.haps$dims) != 2) {
    h5$close_all()
    CacheHaplotypes.err(glue("HDF5 object at '{haps}' does not contain 2-dimensional data."))
  }

  # Check/assign indices to read
  if(!identical(loci.idx, NULL)) {
    if(!transpose) {
      if(!testIntegerish(loci.idx, lower = 1, upper = h5.haps$dims[1], sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, upper = h5.haps$dims[1], sorted = TRUE)}"))
      }
    } else {
      if(!testIntegerish(loci.idx, lower = 1, upper = h5.haps$dims[2], sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, upper = h5.haps$dims[2], sorted = TRUE)}"))
      }
    }
    loci.idx <- loci.idx
  } else {
    if(!transpose) {
      loci.idx <- 1:(h5.haps$dims[1])
    } else {
      loci.idx <- 1:(h5.haps$dims[2])
    }
  }

  if(!identical(hap.idx, NULL)) {
    if(!transpose) {
      if(!testIntegerish(hap.idx, lower = 1, upper = h5.haps$dims[2], sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, upper = h5.haps$dims[2], sorted = TRUE)}"))
      }
    } else {
      if(!testIntegerish(hap.idx, lower = 1, upper = h5.haps$dims[1], sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, upper = h5.haps$dims[1], sorted = TRUE)}"))
      }
    }
    hap.idx <- hap.idx
  } else {
    if(!transpose) {
      hap.idx <- 1:(h5.haps$dims[2])
    } else {
      hap.idx <- 1:(h5.haps$dims[1])
    }
  }

  L <- length(loci.idx)
  N <- length(hap.idx)

  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("N", envir = pkgVars))) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }
  assign("N", as.integer(N), envir = pkgVars)

  # Create a closure for easy access of HDF5 file
  # We'll read in 10MB (raw size) chunks
  make.hdf5.access <- function(h5.haps, loci.idx, hap.idx, transpose, L, N) {
    step.size <- round(10000000/(L/8))
    current.step <- 1
    function() {
      if(current.step > N) {
        return(matrix(nrow = 0, ncol = 0))
      }
      upto <- min(current.step + step.size - 1, N)
      if(!transpose)
        res <- h5.haps[loci.idx,hap.idx[current.step:upto]]
      else
        res <- t(h5.haps[hap.idx[current.step:upto],loci.idx])
      current.step <<- upto + 1
      res
    }
  }

  # Cache it!
  assign("L", CacheHaplotypes_hdf5_2(make.hdf5.access(h5.haps, loci.idx, hap.idx, transpose, L, N), N, L), envir = pkgVars)
  h5$close_all()
  invisible(NULL)
}



CacheHaplotypes.hdf5.rhdf5 <- function(hdf5.file,
                                       loci.idx,
                                       hap.idx,
                                       transpose,
                                       haps) {
  # Check for file and dataset within file
  if(!testFile(hdf5.file, access = "r")) {
    CacheHaplotypes.err(checkFile(hdf5.file, access = "r"))
  }

  # Get dimensions of haplotype data
  hdf5.dim <- integer(2)
  hdf5.dim[1] <- dim(rhdf5::h5read(hdf5.file, haps, index = list(NULL,1)))[1]
  hdf5.dim[2] <- dim(rhdf5::h5read(hdf5.file, haps, index = list(1,NULL)))[2]

  # Check/assign indices to read
  if(!identical(loci.idx, NULL)) {
    if(!transpose) {
      if(!testIntegerish(loci.idx, lower = 1, upper = hdf5.dim[1], sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, upper = hdf5.dim[1], sorted = TRUE)}"))
      }
    } else {
      if(!testIntegerish(loci.idx, lower = 1, upper = hdf5.dim[2], sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument loci.idx ... {checkIntegerish(loci.idx, lower = 1, upper = hdf5.dim[2], sorted = TRUE)}"))
      }
    }
    loci.idx <- loci.idx
  } else {
    if(!transpose) {
      loci.idx <- 1:(hdf5.dim[1])
    } else {
      loci.idx <- 1:(hdf5.dim[2])
    }
  }

  if(!identical(hap.idx, NULL)) {
    if(!transpose) {
      if(!testIntegerish(hap.idx, lower = 1, upper = hdf5.dim[2], sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, upper = hdf5.dim[2], sorted = TRUE)}"))
      }
    } else {
      if(!testIntegerish(hap.idx, lower = 1, upper = hdf5.dim[1], sorted = TRUE)) {
        CacheHaplotypes.err(glue("Error for argument hap.idx ... {checkIntegerish(hap.idx, lower = 1, upper = hdf5.dim[1], sorted = TRUE)}"))
      }
    }
    hap.idx <- hap.idx
  } else {
    if(!transpose) {
      hap.idx <- 1:(hdf5.dim[2])
    } else {
      hap.idx <- 1:(hdf5.dim[1])
    }
  }

  L <- length(loci.idx)
  N <- length(hap.idx)

  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("N", envir = pkgVars))) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }
  assign("N", as.integer(N), envir = pkgVars)

  # Create a closure for easy access of HDF5 file
  # We'll read in 10MB (raw size) chunks
  make.hdf5.access <- function(hdf5.file, loci.idx, hap.idx, transpose, L, N, haps) {
    step.size <- round(10000000/(L/8))
    current.step <- 1
    function() {
      if(current.step > N) {
        return(matrix(nrow = 0, ncol = 0))
      }
      upto <- min(current.step + step.size - 1, N)
      if(!transpose)
        res <- rhdf5::h5read(hdf5.file, haps, index = list(loci.idx, hap.idx[current.step:upto]))
      else
        res <- t(rhdf5::h5read(hdf5.file, haps, index = list(hap.idx[current.step:upto], loci.idx)))
      current.step <<- upto + 1
      res
    }
  }

  # Cache it!
  assign("L", as.integer(CacheHaplotypes_hdf5_2(make.hdf5.access(hdf5.file, loci.idx, hap.idx, transpose, L, N, haps), N, L)), envir = pkgVars)
  invisible(NULL)
}
