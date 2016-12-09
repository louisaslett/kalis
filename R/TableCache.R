# Cache size in GB
CreateForwardTableCache <- function(size = 1, from_recipient = 1, to_recipient = Inf) {
  seqs <- get("seqs", envir = pkgCache)
  if(anyNA(seqs)) {
    stop("No sequences cached ... cannot determine table size until cache is loaded with CacheAllSequences().")
  }
  N <- length(seqs)
  L <- get("seq_size", envir = pkgCache)

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

  cat("Found", N, "sequences in the cache.")
  if((delN*N+2*delN+1)*8/1e9 > size) {
    stop(size, "GB is not big enough for even 1 table.")
  }
  cat("  Constructing table cache of appropriate size ...\n")

  cache <- list()
  i <- 1
  while((length(cache) == 0 || ((object.size(cache)*(length(cache)+1))/length(cache))/1e9 < size) && length(cache)<floor(log2(L))) {
    cache[[i]] <- MakeForwardTable(from_recipient, to_recipient)
    i <- i+1
  }
  cat("Cache constructed, can hold ", length(cache), " tables for recipients ", from_recipient, " ... ", to_recipient, ".  Actual size ≈ ", ceiling(((object.size(cache)*(length(cache)+1))/length(cache))/1e6)/1e3, "GB.\n", sep = "")

  cache
}

ForwardUsingTableCache <- function(fwd, cache, t, Pi, mu, rho, nthreads) {
  l <- sapply(cache, function(x) { x$l })
  if(any(l==t)) {
    CopyForwardTable(fwd, cache[[which(l==t)]])
    return();
  }
  todo <- which(l>t)
  l[todo] <- -1
  from <- max(l)
  from.idx <- which.max(l)

  # Now figure out how to fill in
  if(t-from-1 < length(todo)) {
    # In here, we have more spare slots than there are steps to reach t,
    # so do one step at a time and store into a cache element

    for(i in 1:(t-from-1)) {
      CopyForwardTable(cache[[todo[i]]], cache[[from.idx]])
      Forward(cache[[todo[i]]], from+i, Pi, mu, rho, nthreads)
      from.idx <- todo[i]
    }
    CopyForwardTable(fwd, cache[[from.idx]])
    Forward(fwd, t, Pi, mu, rho, nthreads)
  } else {
    # We have more steps than spare slots, so we need a schedule to fill in the
    # gaps

    b <- exp(log(t-from-1)/length(todo))
    fillin <- unique(floor((1-1/b^{1:length(todo)})*(t-from-1)+from))

    for(i in 1:length(fillin)) {
      CopyForwardTable(cache[[todo[i]]], cache[[from.idx]])
      Forward(cache[[todo[i]]], fillin[i], Pi, mu, rho, nthreads)
      from.idx <- todo[i]
    }
    CopyForwardTable(fwd, cache[[from.idx]])
    Forward(fwd, t, Pi, mu, rho, nthreads)
  }
}
