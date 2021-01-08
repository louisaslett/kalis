##### Gold Master #####
# These are "gold master" implementations of the forward and backward equations
# described in the original kalis paper.  They are "gold master" in the sense
# that they implement the recursions in the simplest (possible least efficient)
# way possible allowing unit testing of the complex C and vector intrinsics code
# against a known correct implementation.  Note, exact equivalence is *not*
# expected due to the vagaries of finite precision floating point arithmetic,
# especially where operations may be reordered in the the high-performance C
# versions.

Theta.GM <- function(pars,locus){
  mu <- ifelse(length(pars$pars$mu)>1, pars$pars$mu[locus], pars$pars$mu[1])

  H <- t(QueryCache(loci.idx = locus))
  if(pars$pars$use.speidel) {
    indH <- outer(H[,1], 1-H[,1], "|")*1; diag(indH) <- 0
  } else {
    indH <- outer(H[,1], H[,1], "==")*1; diag(indH) <- 0
  }

  theta <- indH*(1-2*mu) + mu
}


Forward.GM <- function(fwd, pars, t) {
  res <- MakeForwardTable(pars, fwd$from_recipient, fwd$to_recipient)
  CopyTable(res, fwd)
  res$l = ifelse(fwd$l > L(), 0L, fwd$l)

  while(res$l < t) {
    res$l <- res$l+1L

    if(res$l == 1) {
      res$alpha <- Theta.GM(pars,res$l) * pars$pars$Pi
    } else {
      rho <- pars$pars$rho[res$l-1]
      res$alpha <- Theta.GM(pars,res$l) * (pars$pars$Pi*rho + sweep(res$alpha, 2, (1-rho)/colSums(res$alpha), "*"))
    }
    diag(res$alpha) <- 0 # strictly only needed if scalar Pi is being used
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


  res <- MakeBackwardTable(pars, bck$from_recipient, bck$to_recipient)
  CopyTable(res, bck)
  res$l = ifelse(bck$l > L(), L()+1L, bck$l)

  while(res$l > t) {
    res$l <- res$l-1L

    if(res$l == L()) {
      res$beta <- matrix(1, nrow = N(), ncol = N())
    } else {
      if(!(res$l == (bck$l - 1) & bck$beta.theta)){ # if it's the first backward step off from a table in beta.theta space, beta is already in beta.theta space
                                                    # so this mutation step is already done.
        res$beta <- Theta.GM(pars,res$l+1) * res$beta
      }
      rho <- pars$pars$rho[res$l]
      res$beta <- (1-rho) * sweep(res$beta, 2, colSums(res$beta * pars$pars$Pi), "/") + rho
    }
    diag(res$beta) <- 0
  }

  res$beta.g <- colSums(Theta.GM(pars,res$l) * res$beta * pars$pars$Pi) # FIX ME! beta.g is not right here

  if(beta.theta){
    res$beta <- Theta.GM(pars,res$l) * res$beta
    res$beta.theta <- TRUE
  }

  CopyTable(bck, res)

  invisible(NULL)
}
