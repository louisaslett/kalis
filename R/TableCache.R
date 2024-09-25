#' Create cache for forward tables
#'
#' Create an in-memory cache for forward tables to improve efficiency when
#' iterating in reverse along the haplotype sequences.
#'
#' If the objective is to run the Li and Stephens hidden Markov model both
#' forwards and backwards to the same locus (and to do so for every possible
#' locus), then considerable efficiency can be achieved by first performing a
#' full scan forwards, filling a geometrically spaced cache whilst doing so.
#' Then, by working backwards, the backward propagation moves one locus at a
#' time and the forward propagation can move backwards by moving forward from a
#' recently cached local table.
#'
#' Memory for a cache can be allocated using this function and should then be
#' filled with [FillTableCache()].
#' To use the cache, then instead of using the [Forward()] function,
#' use [ForwardUsingTableCache()].
#'
#' @param pars a `kalisParameters` object, as returned by [Parameters()].
#' @param size the maximum amount of RAM (in GB) to devote to this cache.
#' @param from_recipient first recipient haplotype if creating a partial forward
#'   table cache.  By default includes from the first recipient haplotype.
#' @param to_recipient last recipient haplotype if creating a partial forward
#'   table cache.  By default includes to the last recipient haplotype.
#'
#' @return
#'   A list of forward tables representing a cache and ready to be filled is
#'   returned.
#'
#' @seealso
#'   [MakeForwardTable()] to make a forward table;
#'   [FillTableCache()] to fill a cache;
#'   [ForwardUsingTableCache()] to use a cache;
#'   [Forward()] for forward function without using a cache.
#'
#' @examples
#' \dontrun{
#' # This code assumes you have already:
#' #  i) cached the haplotypes using CacheHaplotypes function
#' #  ii) setup parameters in a variable called pars
#' #  iii) set the number of loci in a variable called L
#'
#' # Allocate up to 10GB to a cache, with parameters already setup in pars ...
#' cache <- CreateForwardTableCache(pars, 10)
#' # ... and fill it
#' FillTableCache(cache, pars, nthreads = 8)
#'
#' # Create forward and backward tables
#' fwd <- MakeForwardTable(pars)
#' bck <- MakeBackwardTable(pars)
#'
#' # Then reach every locus faster by iterating backwards, using the cache to
#' # move the forward table into position faster
#' for(l in L:1) {
#'   Backward(bck, pars, l, nthreads = 8)
#'   ForwardUsingTableCache(fwd, pars, cache, l, nthreads = 8)
#'   # Do whatever work is required at
#'   # every locus here using fwd and bck
#' }
#' }
#'
#' @export
CreateForwardTableCache <- function(pars, size = 1, from_recipient = 1, to_recipient = Inf, max.tables = 0) {
  if(!("kalisParameters" %in% class(pars))) {
    stop("The pars argument is not a valid parameters object.")
  }

  N <- get("N", envir = pkgVars)
  L <- get("L", envir = pkgVars)

  if(anyNA(N)) {
    stop("No haplotypes cached ... cannot determine table size until cache is loaded with CacheAllHaplotypes().")
  }

  if(from_recipient>to_recipient) {
    stop("from_recipient must be smaller than to_recipient.")
  }
  if(from_recipient < 1) {
    from_recipient <- 1
  }
  if(to_recipient > N) {
    to_recipient <- N
  }
  delN <- to_recipient-from_recipient+1
  if(!is.vector(max.tables) || !is.numeric(max.tables) || length(max.tables) != 1 || max.tables < 0) {
    stop("max.tables must be a positive scalar.")
  }

  cat("Found", N, "haplotypes in the cache.")
  if((delN*N+2*delN+1)*8/1e9 > size) {
    stop(size, "GB is not big enough for even 1 table.")
  }
  cat("  Constructing table cache of appropriate size ...\n")

  if(max.tables == 0) {
    max.tables <- floor(log2(L))
  }
  cache <- list()
  i <- 1
  while((length(cache) == 0 || ((utils::object.size(cache)*(length(cache)+1))/length(cache))/1e9 < size) && length(cache)<max.tables) {
    cache[[i]] <- MakeForwardTable(pars, from_recipient, to_recipient)
    i <- i+1
  }
  cat("Cache constructed, can hold ", length(cache), " tables for recipients ", from_recipient, " ... ", to_recipient, ".  Actual size approx ", ceiling(utils::object.size(cache)/1e6)/1e3, "GB.\n", sep = "")

  class(cache) <- c("kalisCheckpointTable", class(cache))
  cache
}



#' Fill a forward table cache
#'
#' An in-memory cache for forward tables can be filled using this function, for
#' either the whole sequence length or some sub-sequence.
#'
#' If the objective is to run the Li and Stephens hidden Markov model both
#' forwards and backwards to the same locus (and to do so for every possible
#' locus), then considerable efficiency can be achieved by first performing a
#' full scan forwards, filling a geometrically spaced cache whilst doing so.
#' Then, by working backwards, the backward propagation moves one locus at a
#' time and the forward propagation can move backwards by moving forward from a
#' recently cached local table.
#'
#' Memory for a cache can be allocated using
#' [CreateForwardTableCache()] and should then be filled with this
#' function.
#' To use the cache, then instead of using the [Forward()] function,
#' use [ForwardUsingTableCache()].
#'
#' @param cache a cache of forward tables as generated by
#'   [CreateForwardTableCache()]
#' @param pars a `kalisParameters` object, as returned by [Parameters()].
#' @param from the first locus which the geometrically spaced cache should be
#'   built from.
#'   By default the whole sequence length will be cached so this defaults to 1.
#' @param to the last locus up to which the geometrically spaced cache should be
#'   built.
#'   By default the whole sequence length will be cached so this defaults to
#'   `Inf`.
#' @param nthreads the number of CPU cores to use.
#'   By default no parallelism is used.
#'
#' @return
#' There is nothing returned.
#' For performance reasons, `cache` is updated in-place.
#'
#' @seealso
#'   [MakeForwardTable()] to make a forward table;
#'   [CreateForwardTableCache()] to generate a cache;
#'   [ForwardUsingTableCache()] to use a cache;
#'   [Forward()] for forward function without using a cache.
#'
#' @examples
#' \dontrun{
#' # This code assumes you have already:
#' #  i) cached the haplotypes using CacheHaplotypes function
#' #  ii) setup parameters in a variable called pars
#' #  iii) set the number of loci in a variable called L
#'
#' # Allocate up to 10GB to a cache, with parameters already setup in pars ...
#' cache <- CreateForwardTableCache(pars, 10)
#' # ... and fill it
#' FillTableCache(cache, pars, nthreads = 8)
#'
#' # Create forward and backward tables
#' fwd <- MakeForwardTable(pars)
#' bck <- MakeBackwardTable(pars)
#'
#' # Then reach every locus faster by iterating backwards, using the cache to
#' # move the forward table into position faster
#' for(l in L:1) {
#'   Backward(bck, pars, l, nthreads = 8)
#'   ForwardUsingTableCache(fwd, pars, cache, l, nthreads = 8)
#'   # Do whatever work is required at
#'   # every locus here using fwd and bck
#' }
#' }
#'
#' @export
FillTableCache <- function(cache,
                           pars,
                           nthreads = min(parallel::detectCores(logical = FALSE), cache[[1]]$to_recipient-cache[[1]]$from_recipient+1)) {
  if(!("kalisCheckpointTable" %in% class(cache))) {
    stop("The cache argument is not a valid forward table cache.")
  }
  if(!("kalisParameters" %in% class(pars))) {
    stop("The pars argument is not a valid parameters object.")
  }
  if(cache[[1]]$pars.sha256 != pars$sha256) {
    stop("The forward table provided was created with different parameter values (SHA-256 mismatch).")
  }

  L <- get("L", envir = pkgVars)

  from <- cache[[1]]$from_recipient
  to <- cache[[1]]$from_recipient

  if(from < 1 || from > L) {
    stop("from argument is invalid.")
  }
  if(is.infinite(to)) {
    to <- L
  }
  if(to < 1 || to > L) {
    stop("to argument is invalid.")
  }

  pos <- 0.5
  if(from > 1) {
    pos <- 1
  }

  for(i in 1:length(cache)) {
    t <- floor((1.0-pos)*(to-from))+from
    if(t == cache[[i]]$l) {
      ResetTable(cache[[i]])
      break
    }
    pos <- pos*0.5

    cat(glue("Computing cache entry {i} up to locus {t} for recipients {cache[[i]]$from_recipient} to {cache[[i]]$to_recipient} from "))
    if(cache[[i]]$l < 1) {
      cat("start\n")
    } else {
      cat(glue("locus {cache[[i]]$l}"), "\n")
    }

    Forward(cache[[i]], pars, t, nthreads)

    if(i < length(cache)) {
      CopyTable(cache[[i+1]], cache[[i]])
    }
  }
}



#' Use a forward table cache to propagate
#'
#' An in-memory cache for forward tables, which has already been filled, can be
#' used to move more quickly to a specified locus.
#'
#' If the objective is to run the Li and Stephens hidden Markov model both
#' forwards and backwards to the same locus (and to do so for every possible
#' locus), then considerable efficiency can be achieved by first performing a
#' full scan forwards, filling a geometrically spaced cache whilst doing so.
#' Then, by working backwards, the backward propagation moves one locus at a
#' time and the forward propagation can move backwards by moving forward from a
#' recently cached local table.
#'
#' Memory for a cache can be allocated using
#' [CreateForwardTableCache()] and should then be filled with
#' [FillTableCache()].
#' To use the cache, then instead of using the [Forward()] function,
#' use this function.
#'
#' Note that the `cache` which is passed to this function will be
#' dynamically updated based on the locus requested: the assumption is that
#' the cache is used to propagate in reverse so any cache entries for a locus
#' position past `t` are taken to be no longer needed and that space will
#' redeployed to more densely fill the cache with earlier locus positions.
#'
#' @param fwd a `kalisForwardTable` object, as returned by
#'   [MakeForwardTable()].
#' @param pars a `kalisParameters` object, as returned by [Parameters()].
#' @param cache a cache of forward tables as generated by
#'   [CreateForwardTableCache()] and filled using
#'   [FillTableCache()].
#' @param t a locus position to move the forward table to, starting the forward
#'   propagation from whatever table in the `cache` variable is immediately
#'   before locus `t`.
#'   By default, it simply advances to the previous locus (which is the natural
#'   direction to move when using the cache).
#' @param nthreads the number of CPU cores to use.
#'   By default no parallelism is used.
#'
#' @return
#'   There is nothing returned.
#'   For performance reasons, `fwd` is updated in-place.
#'
#' @seealso
#'   [MakeForwardTable()] to make a forward table;
#'   [CreateForwardTableCache()] to generate a cache;
#'   [FillTableCache()] to fill a cache;
#'   [Forward()] for forward function without using a cache.
#'
#' @examples
#' \dontrun{
#' # This code assumes you have already:
#' #  i) cached the haplotypes using CacheHaplotypes function
#' #  ii) setup parameters in a variable called pars
#' #  iii) set the number of loci in a variable called L
#'
#' # Allocate up to 10GB to a cache, with parameters already setup in pars ...
#' cache <- CreateForwardTableCache(pars, 10)
#' # ... and fill it
#' FillTableCache(cache, pars, nthreads = 8)
#'
#' # Create forward and backward tables
#' fwd <- MakeForwardTable(pars)
#' bck <- MakeBackwardTable(pars)
#'
#' # Then reach every locus faster by iterating backwards, using the cache to
#' # move the forward table into position faster
#' for(l in L:1) {
#'   Backward(bck, pars, l, nthreads = 8)
#'   ForwardUsingTableCache(fwd, pars, cache, l, nthreads = 8)
#'   # Do whatever work is required at
#'   # every locus here using fwd and bck
#' }
#' }
#'
#' @export
ForwardUsingTableCache <- function(fwd, pars, cache, t = fwd$l-1,
                                   nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient-fwd$from_recipient+1)) {
  if(!("kalisForwardTable" %in% class(fwd))) {
    stop("The fwd argument is not a valid forward table.")
  }
  if(!("kalisParameters" %in% class(pars))) {
    stop("The pars argument is not a valid parameters object.")
  }
  if(fwd$pars.sha256 != pars$sha256) {
    stop("The forward table provided was created with different parameter values (SHA-256 mismatch).")
  }
  if(!("kalisCheckpointTable" %in% class(cache))) {
    stop("The cache argument is not a valid forward table cache.")
  }
  if(cache[[1]]$pars.sha256 != pars$sha256) {
    stop("The forward table cache provided was created with different parameter values (SHA-256 mismatch).")
  }

  if(!all(sapply(cache,function(x){x$from_recipient}) == fwd$from_recipient)){
    stop("from_recipient is not consistent across fwd and the tables in the cache.")
  }

  if(!all(sapply(cache,function(x){x$to_recipient}) == fwd$to_recipient)){
    stop("to_recipient is not consistent across fwd and the tables in the cache.")
  }


  l <- sapply(cache, function(x) { x$l })
  if(any(l==t)) {
    CopyTable(fwd, cache[[which(l==t)]])
    return()
  }
  todo <- which(l>t)
  l[todo] <- -1
  # Is the max -1?  Then we've passed the half way mark.  Fill up the first slot
  # and march on
  if(max(l) == -1) {
    ResetTable(cache[[1]]) # Have to do this in C++ as must in place modify
    Forward(cache[[1]], pars, 1, nthreads)
    l[1] <- 1
    todo <- which(l==-1)
  }
  from <- max(l)
  from.idx <- which.max(l)

  # First, check if there are any spare slots -- we might just be accessing after
  # the last checkpoint already.  If so, run forward and return right away
  if(length(todo)==0) {
    CopyTable(fwd, cache[[from.idx]])
    Forward(fwd, pars, t, nthreads)
    return()
  }

  # If we want just one step after the jumping off checkpoint, then we just wind
  # forward to it right away though I have spare checkpoint slots
  if(t == from+1) {
    CopyTable(fwd, cache[[from.idx]])
    Forward(fwd, pars, t, nthreads)
    return()
  }

  # Now figure out how to fill in
  if(t-from <= length(todo)) {
    # In here, we have more spare slots than there are steps to reach t,
    # so do one step at a time and store into a cache element

    # NB t-from >= 2 due to if statement above
    for(i in 1:(t-from-1)) {
      CopyTable(cache[[todo[i]]], cache[[from.idx]])
      Forward(cache[[todo[i]]], pars, from+i, nthreads)
      from.idx <- todo[i]
    }
    CopyTable(fwd, cache[[from.idx]])
    Forward(fwd, pars, t, nthreads)
  } else {
    # We have more steps than spare slots, so we need a schedule to fill in the
    # gaps

    # NB t-from-1 >= 1 due to if statement above
    geom.spacing <- ceiling( (1-0.5^{1:length(todo)})*(t-from-1)+from )
    fillin <- rep(NA, length(geom.spacing))

    fillin[length(geom.spacing)] <- geom.spacing[length(geom.spacing)]
    if(length(geom.spacing)>1) {
      for(i in (length(geom.spacing)-1):1) {
        fillin[i] <- min(geom.spacing[i], fillin[i+1]-1)
      }
    }
    fillin <- unique(fillin[fillin>from & fillin<t])

    for(i in 1:length(fillin)) {
      CopyTable(cache[[todo[i]]], cache[[from.idx]])
      Forward(cache[[todo[i]]], pars, fillin[i], nthreads)
      from.idx <- todo[i]
    }
    CopyTable(fwd, cache[[from.idx]])
    Forward(fwd, pars, t, nthreads)
  }
}



#' @export
print.kalisCheckpointTable <- function(x, ...) {
  if(!("kalisCheckpointTable" %in% class(x))) stop("Not a kalisCheckpointTable object")
  cat("Checkpoint Table object containing", length(x), "checkpoints.\n")
  cat("  Loci of checkpoints:\n")
  cat("   ", sapply(x, function(x) { x$l }), "\n")
  cat("  Memory consumed: ", ceiling(utils::object.size(x)/1e6)/1e3, "GB.\n")
}