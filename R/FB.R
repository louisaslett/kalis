Forward <- function(fwd, t, Pi, mu, rho, nthreads) {
  L <- get("seq_size", envir = pkgCache)
  N <- length(get("seqs", envir = pkgCache))
  if(fwd$l > t) {
    stop("The forward table provided is for locus position ", l, " which is already past requested locus ", t)
  }
  if(nrow(fwd$alpha) != N || ncol(fwd$alpha) != fwd$to_recipient-fwd$from_recipient+1) {
    stop("Forward table is of the wrong dimensions for this problem.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.atomic(Pi) && length(Pi) == 1L && !is.character(Pi) && Im(Pi)==0)) {
    stop("If Pi is not a full matrix of background copying probabilities then it must be a single scalar for uniform background copying probability.")
  }
  if(!is.vector(mu) || length(mu) != L) {
    stop("mu is of the wrong type/length.")
  }
  if(!is.vector(rho) || length(rho) != L) {
    stop("rho is of the wrong type/length.")
  }

  if(is.matrix(Pi)) {
    Forward_densePi_cpp(fwd, t, Pi, mu, rho, nthreads)
  } else {
    Forward_scalarPi_cpp(fwd, t, Pi, mu, rho, nthreads)
  }
}

Backward <- function(bck, t, Pi, mu, rho, nthreads) {
  L <- get("seq_size", envir = pkgCache)
  N <- length(get("seqs", envir = pkgCache))
  if(bck$l < t) {
    stop("The backward table provided is for locus position ", l, " which is already before requested locus ", t)
  }
  if(nrow(bck$beta) != N || ncol(bck$beta) != bck$to_recipient-bck$from_recipient+1) {
    stop("Forward table is of the wrong dimensions for this problem.")
  }
  if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
    stop("Pi is of the wrong dimensions for this problem.")
  }
  if(!is.matrix(Pi) && !(is.atomic(Pi) && length(Pi) == 1L && !is.character(Pi) && Im(Pi)==0)) {
    stop("If Pi is not a full matrix of background copying probabilities then it must be a single scalar for uniform background copying probability.")
  }
  if(!is.vector(mu) || length(mu) != L) {
    stop("mu is of the wrong type/length.")
  }
  if(!is.vector(rho) || length(rho) != L) {
    stop("rho is of the wrong type/length.")
  }

  if(is.matrix(Pi)) {
    Backward_densePi_cpp(bck, t, Pi, mu, rho, nthreads)
  } else {
    Backward_scalarPi_cpp(bck, t, Pi, mu, rho, nthreads)
  }
}





