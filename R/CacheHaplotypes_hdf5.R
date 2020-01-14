CacheHaplotypes.hdf5 <- function(hdf5.file,
                                 transpose = FALSE,
                                 haps = "/haps",
                                 hap.ids = NULL,
                                 loci.ids = NULL,
                                 hdf5.pkg = c("hdf5r", "rhdf5")) {
  for(pkg in hdf5.pkg) {
    if(!(pkg %in% c("hdf5r", "rhdf5"))) {
      stop(glue("The package {pkg} is not supported for reading from HDF5 files."))
    }
    if(length(find.package(pkg, quiet = TRUE)) > 0) {
      if(pkg == "hdf5r") {
        return(CacheHaplotypes.hdf5.hdf5r(hdf5.file, transpose, haps, hap.ids, loci.ids))
      }
      if(pkg == "rhdf5") {
        return(CacheHaplotypes.hdf5.rhdf5(hdf5.file, transpose, haps, hap.ids, loci.ids))
      }
    }
  }
}



CacheHaplotypes.hdf5.hdf5r <- function(hdf5.file,
                                       transpose,
                                       haps,
                                       hap.ids,
                                       loci.ids) {
  # Check for file and dataset within file
  if(!file.exists(hdf5.file)) {
    CacheHaplotypes.err("Cannot find HDF5 file.")
  }

  # Open file and objects
  h5 <- hdf5r::H5File$new(hdf5.file, "r")
  if(!h5$exists(haps)) {
    h5$close_all()
    CacheHaplotypes.err(glue("HDF5 file exists but does not contain an object at '{haps}' for the haplotype data."))
  }

  # haplotypes
  h5.haps <- h5$open(haps)
  if(length(h5.haps$dims) != 2) {
    h5$close_all()
    CacheHaplotypes.err(glue("HDF5 object at '{haps}' does not contain 2-dimensional data."))
  }
  if(!transpose) {
    N <- h5.haps$dims[2]
    L <- h5.haps$dims[1]
  } else {
    N <- h5.haps$dims[1]
    L <- h5.haps$dims[2]
  }

  # haplotype ids
  if(!is.null(hap.ids)) {
    if(!h5$exists(hap.ids)) {
      h5$close_all()
      CacheHaplotypes.err(glue("HDF5 file exists but does not contain an object at '{hap.ids}' for the haplotype IDs."))
    }
    h5.hap.ids <- h5$open(hap.ids)
    if(length(h5.hap.ids$dims) != 1) {
      h5$close_all()
      CacheHaplotypes.err(glue("HDF5 object at '{hap.ids}' does not contain 1-dimensional ID information."))
    }
    if(h5.hap.ids$dims[1] != N) {
      h5$close_all()
      CacheHaplotypes.err(glue("HDF5 object at '{hap.ids}' contains {h5.hap.ids$dims[1]} haplotype IDs, but {haps} dataset contains {N} haplotypes (is transpose set correctly?)."))
    }
    hap.ids.tmp <- as.character(h5.hap.ids[])
    if(length(unique(hap.ids.tmp)) != length(hap.ids.tmp)) {
      h5$close_all()
      CacheHaplotypes.err(glue("HDF5 object at '{hap.ids}' contains {h5.hap.ids$dims[1]} haplotype IDs, but they are not all unique."))
    }
  } else {
    hap.ids.tmp <- as.integer(1:N)
  }

  # loci ids
  if(!is.null(loci.ids)) {
    if(!h5$exists(loci.ids)) {
      h5$close_all()
      CacheHaplotypes.err(glue("HDF5 file exists but does not contain an object at '{loci.ids}' for the locus IDs."))
    }
    h5.loci.ids <- h5$open(loci.ids)
    if(length(h5.loci.ids$dims) != 1) {
      h5$close_all()
      CacheHaplotypes.err(glue("HDF5 object at '{loci.ids}' does not contain 1-dimensional ID information."))
    }
    if(h5.loci.ids$dims[1] != L) {
      h5$close_all()
      CacheHaplotypes.err(glue("HDF5 object at '{loci.ids}' contains {h5.loci.ids$dims[1]} locus IDs, but {haps} dataset contains {L} loci (is transpose set correctly?)."))
    }
    loci.ids.tmp <- as.character(h5.loci.ids[])
    if(length(unique(loci.ids.tmp)) != length(loci.ids.tmp)) {
      h5$close_all()
      CacheHaplotypes.err(glue("HDF5 object at '{loci.ids}' contains {h5.loci.ids$dims[1]} locus IDs, but they are not all unique."))
    }
  } else {
    loci.ids.tmp <- as.integer(1:L)
  }

  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("N", envir = pkgVars))) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }
  assign("N", as.integer(N), envir = pkgVars)
  assign("hap.ids", hap.ids.tmp, envir = pkgVars); rm(hap.ids.tmp);
  assign("loci.ids", loci.ids.tmp, envir = pkgVars); rm(loci.ids.tmp);

  # Create a closure for easy access of HDF5 file
  # We'll read in 10MB (raw size) chunks
  make.hdf5.access <- function(h5.haps, N, L) {
    step.size <- round(10000000/(L/8))
    current.step <- 1
    function() {
      if(current.step > N) {
        return(matrix(nrow = 0, ncol = 0))
      }
      upto <- min(current.step + step.size - 1, N)
      if(!transpose)
        res <- h5.haps[,current.step:upto]
      else
        res <- t(h5.haps[current.step:upto,])
      current.step <<- upto + 1
      res
    }
  }

  # Cache it!
  assign("L", CacheHaplotypes_hdf5_2(make.hdf5.access(h5.haps, N, L), N, L), envir = pkgVars)
  h5$close_all()
  invisible(NULL)
}



CacheHaplotypes.hdf5.rhdf5 <- function(hdf5.file,
                                       transpose,
                                       haps,
                                       hap.ids,
                                       loci.ids) {
  # Check for file and dataset within file
  if(!file.exists(hdf5.file)) {
    CacheHaplotypes.err("Cannot find HDF5 file.")
  }

  # Get dimensions of haplotype data
  hdf5.dim <- integer(2)
  hdf5.dim[1] <- dim(rhdf5::h5read(hdf5.file, haps, index = list(NULL,1)))[1]
  hdf5.dim[2] <- dim(rhdf5::h5read(hdf5.file, haps, index = list(1,NULL)))[2]
  if(!transpose) {
    N <- hdf5.dim[2]
    L <- hdf5.dim[1]
  } else {
    N <- hdf5.dim[1]
    L <- hdf5.dim[2]
  }

  # Haplotype and locus IDs
  if(is.null(hap.ids)) {
    hap.ids.tmp <- as.integer(1:N)
  } else {
    hap.ids.tmp <- as.character(rhdf5::h5read(hdf5.file, hap.ids))
  }
  if(is.null(loci.ids)) {
    loci.ids.tmp <- as.integer(1:L)
  } else {
    loci.ids.tmp <- as.character(rhdf5::h5read(hdf5.file, loci.ids))
  }
  if(length(hap.ids.tmp) != N) {
    CacheHaplotypes.err(glue("There are {length(hap.ids.tmp)} haplotype IDs provided at '{hap.ids}' in the HDF5 file, but there are {N} haplotypes provided in '{haps}' (do you need to change transpose?)"))
  }
  if(length(unique(hap.ids.tmp)) != length(hap.ids.tmp)) {
    CacheHaplotypes.err(glue("HDF5 object at '{hap.ids}' contains {length(hap.ids.tmp)} haplotype IDs, but they are not all unique."))
  }
  if(length(loci.ids.tmp) != L) {
    CacheHaplotypes.err(glue("There are {length(loci.ids.tmp)} locus IDs provided at '{loci.ids}' in the HDF5 file, but there are {L} loci provided in '{haps}' (do you need to change transpose?)"))
  }
  if(length(unique(loci.ids.tmp)) != length(loci.ids.tmp)) {
    CacheHaplotypes.err(glue("HDF5 object at '{loci.ids}' contains {length(loci.ids.tmp)} locus IDs, but they are not all unique."))
  }

  # Make sure we don't double up cache content if there is already stuff cached
  if(!is.na(get("N", envir = pkgVars))) {
    warning("haplotypes already cached ... overwriting existing cache.")
    ClearHaplotypeCache()
  }
  assign("N", as.integer(N), envir = pkgVars)
  assign("hap.ids", hap.ids.tmp, envir = pkgVars); rm(hap.ids.tmp);
  assign("loci.ids", loci.ids.tmp, envir = pkgVars); rm(loci.ids.tmp);

  # Create a closure for easy access of HDF5 file
  # We'll read in 10MB (raw size) chunks
  make.hdf5.access <- function(hdf5.file, N, L, haps) {
    step.size <- round(10000000/(L/8))
    current.step <- 1
    function() {
      if(current.step > N) {
        return(matrix(nrow = 0, ncol = 0))
      }
      upto <- min(current.step + step.size - 1, N)
      if(!transpose)
        res <- rhdf5::h5read(hdf5.file, haps, index = list(NULL, current.step:upto))
      else
        res <- t(rhdf5::h5read(hdf5.file, haps, index = list(current.step:upto, NULL)))
      current.step <<- upto + 1
      res
    }
  }

  # Cache it!
  assign("L", as.integer(CacheHaplotypes_hdf5_2(make.hdf5.access(hdf5.file, N, L, haps), N, L)), envir = pkgVars)
  invisible(NULL)
}
