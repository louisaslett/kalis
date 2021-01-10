##### Gold Master #####
# These are "gold master" implementations of the forward and backward equations
# described in the original kalis paper.  They are "gold master" in the sense
# that they implement the recursions in the simplest (possible least efficient)
# way possible allowing unit testing of the complex C and vector intrinsics code
# against a known correct implementation.  Note, exact equivalence is *not*
# expected due to the vagaries of finite precision floating point arithmetic,
# especially where operations may be reordered in the the high-performance C
# versions.

Theta.GM <- function(pars,locus,from,to){
  mu <- ifelse(length(pars$pars$mu)>1, pars$pars$mu[locus], pars$pars$mu[1])

  H <- c(QueryCache(loci.idx = locus))
  if(pars$pars$use.speidel) {
    indH <- outer(H, 1-H[from:to], "|")*1; diag(indH[from:to,]) <- 0
  } else {
    indH <- outer(H, H[from:to], "==")*1; diag(indH[from:to,]) <- 0
  }

  theta <- indH*(1-2*mu) + mu
}


Forward.GM <- function(fwd, pars, t) {
  res <- MakeForwardTable(pars, fwd$from_recipient, fwd$to_recipient)
  CopyTable(res, fwd)
  res$l = ifelse(fwd$l > L(), 0L, fwd$l)

  if(is.matrix(pars$pars$Pi)) {
    Pi <- pars$pars$Pi[,fwd$from_recipient:fwd$to_recipient]
  } else {
    Pi <- pars$pars$Pi
  }

  while(res$l < t) {
    res$l <- res$l+1L

    if(res$l == 1) {
      res$alpha <- Theta.GM(pars,res$l, fwd$from_recipient, fwd$to_recipient) * Pi
    } else {
      rho <- pars$pars$rho[res$l-1]
      res$alpha <- Theta.GM(pars,res$l, fwd$from_recipient, fwd$to_recipient) * (Pi*rho + sweep(res$alpha, 2, (1-rho)/colSums(res$alpha), "*"))
    }
    diag(res$alpha[fwd$from_recipient:fwd$to_recipient,]) <- 0 # strictly only needed if scalar Pi is being used
  }

  res$alpha.f <- colSums(res$alpha)
  CopyTable(fwd, res)

  invisible(NULL)
}



Backward.GM <- function(bck, pars, t, beta.theta = FALSE) {

  if(bck$l == t & bck$beta.theta == beta.theta){ # Nothing to be done.
    return(invisible(NULL))
  }

  if(bck$l == t & bck$beta.theta & !beta.theta){
    warning("Cannot move from beta.theta space back to beta space while staying at the same locus.")
    return(invisible(NULL))
  }

  if(is.matrix(pars$pars$Pi)) {
    Pi <- pars$pars$Pi[,bck$from_recipient:bck$to_recipient]
  } else {
    Pi <- pars$pars$Pi
  }

  res <- MakeBackwardTable(pars, bck$from_recipient, bck$to_recipient)
  CopyTable(res, bck)
  res$l = ifelse(bck$l > L(), L()+1L, bck$l)

  while(res$l > t) {
    res$l <- res$l-1L

    if(res$l == L()) {
      res$beta <- matrix(1, nrow = N(), ncol = bck$to_recipient - bck$from_recipient + 1)
    } else {
      if(!(res$l == (bck$l - 1) & bck$beta.theta)){ # if it's the first backward step off from a table in beta.theta space, beta is already in beta.theta space
        # so this mutation step is already done.
        res$beta <- Theta.GM(pars,res$l+1, bck$from_recipient, bck$to_recipient) * res$beta
      }
      rho <- pars$pars$rho[res$l]
      res$beta <- (1-rho) * sweep(res$beta, 2, colSums(res$beta * Pi), "/") + rho
    }
    diag(res$beta[bck$from_recipient:bck$to_recipient,]) <- 0
  }

  res$beta.g <- colSums(Theta.GM(pars,res$l, bck$from_recipient, bck$to_recipient) * res$beta * Pi) # FIX ME! beta.g is not right here

  if(beta.theta){
    res$beta <- Theta.GM(pars,res$l, bck$from_recipient, bck$to_recipient) * res$beta
    res$beta.theta <- TRUE
  }

  CopyTable(bck, res)

  invisible(NULL)
}



PostProbs.GM <- function(fwd, bck, unif.on.underflow = FALSE, beta.theta.opts = NULL) {

  N <- nrow(fwd$alpha)

  rho.list <- input_checks_for_probs_and_dist_mat(fwd,bck,beta.theta.opts)

  if(!is.null(beta.theta.opts) && !is.null(beta.theta.opts$pars)){
    Pi <- beta.theta.opts$pars$pars$Pi
  } else {
    Pi <- 1/(N-1)
  }

  if(bck$beta.theta){
    # propagate fwd and bck as done in our standard Forward/Backward recursions
    probs <- (sweep(fwd$alpha,2,fwd$alpha.f,"/") * (1-rho.list$rho.fwd) + rho.list$rho.fwd * Pi) *
      (sweep(bck$beta,2,bck$beta.g,"/") * (1-rho.list$rho.bck) + rho.list$rho.bck)

    probs[!is.finite(probs)] <- 0 # covers all cases of the input being NaN or Inf, or dividing any finite numbers or 0 by 0.
    diag(probs[fwd$from_recipient:fwd$to_recipient,]) <- 0
  } else {
    probs <- fwd$alpha * bck$beta
  }


  z0 <- colSums(probs)
  probs <- sweep(probs,2,z0,"/")

  if(unif.on.underflow){
    probs[!is.finite(probs)] <- 1/(N-1) # covers all cases of the input being NaN or Inf, or dividing any finite numbers or 0 by 0.
    diag(probs[fwd$from_recipient:fwd$to_recipient,]) <- 0
  } else {
    probs[!is.finite(probs)] <- 0
  }

  probs
}


DistMat.GM <- function(fwd, bck, type = "raw", beta.theta.opts = NULL) {

  N <- nrow(fwd$alpha)

  rho.list <- input_checks_for_probs_and_dist_mat(fwd,bck,beta.theta.opts)

  if(!is.null(beta.theta.opts) && !is.null(beta.theta.opts$pars)){
    Pi <- beta.theta.opts$pars$pars$Pi
  } else {
    Pi <- 1/(N-1)
  }


  if(bck$beta.theta){
    # propagate fwd and bck as done in our standard Forward/Backward recursions
    probs <- (sweep(fwd$alpha,2,fwd$alpha.f,"/") * (1-rho.list$rho.fwd) + rho.list$rho.fwd * Pi) *
      (sweep(bck$beta,2,bck$beta.g,"/") * (1-rho.list$rho.bck) + rho.list$rho.bck)

    probs[!is.finite(probs)] <- 0 # covers all cases of the input being NaN or Inf, or dividing any finite numbers or 0 by 0.
    diag(probs[fwd$from_recipient:fwd$to_recipient,]) <- 0
  } else {
    probs <- fwd$alpha * bck$beta
  }

  if(type == "raw" | type == "std"){
    z0 <- colSums(probs)
  } else {
    z0 <- c(apply(probs,2,max))
  }

  probs <- sweep(probs,2,z0,"/")
  dists <- -log(probs)

  dists[!is.finite(dists)] <- 744.4400719213812180897 # covers all cases of the input being NaN or Inf, or dividing any finite numbers or 0 by 0.
  diag(dists[fwd$from_recipient:fwd$to_recipient,]) <- 0

  if(type == "raw" | type == "minus.min"){ return(dists) }

  diag(dists[fwd$from_recipient:fwd$to_recipient,]) <- NA

  dm <- colMeans(dists,na.rm = T)
  ds <- apply(dists,2,sd,na.rm = TRUE)

  dists <- sweep(dists,2,dm,"-")
  dists <- sweep(dists,2,ds,"/")
  dists[!is.finite(dists)] <- 0 # covers dividing any finite numbers or 0 by 0.
  diag(dists[fwd$from_recipient:fwd$to_recipient,]) <- 0

  dists
}

