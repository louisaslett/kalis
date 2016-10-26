##### Gold Master Versions #####
# Gold master exact table computation, matrix-wise, log-space
Exact_ComputeTable_GMmatLS_R <- function(l, Pi, mu, rho) {
  # Load/initialise variables we require
  L <- get("seq_size", envir = pkgCache)

  seqs <- get("seqs", envir = pkgCache)
  N <- length(seqs)

  H <- matrix(0, nrow = N, ncol = L)
  for(i in 1:N) {
    H[i,] <- QueryCache(seqs[i])
  }

  # Sanity checks
  if(!is.matrix(Pi)) {
    stop("Haplotype copy probabilities (Pi) should be a matrix.")
  }
  if(!all(dim(Pi) == N)) {
    stop("Haplotype copy probability matrix (Pi) should be N x N.")
  }
  if(!is.vector(mu)) {
    stop("Mutation rate/penalty (mu) should be a vector.")
  }
  if(length(mu) != L) {
    stop("Mutation rate/penalty (mu) should be length L.")
  }
  if(!is.vector(rho)) {
    stop("Scaled expected number of recombination events (rho) should be a vector.")
  }
  if(length(rho) != L) {
    stop("Scaled expected number of recombination events (rho) should be length L.")
  }

  # Compute table from only GM functions
  alphasum <- Exact_Forward_GMmatLS_R(L, L, N, H, Pi, mu, rho)$alphasum

  res.fwd <- Exact_Forward_GMmatLS_R(l, L, N, H, Pi, mu, rho)$alpha
  res.rev <- Exact_Backward_GMmatLS_R(l, L, N, H, Pi, mu, rho)

  list(fwd=res.fwd, rev=res.rev, alphasum=alphasum)
}
# Gold master backward algorithm, no sanity checks, matrix-wise, log-space
Exact_Backward_GMmatLS_R <- function(t, L, N, H, Pi, mu, rho) {
  # Setup Pi & rho
  diag(Pi) <- 0

  l <- L

  indH <- outer(H[,l], H[,l], "==") * 1; diag(indH) <- 0
  theta <- indH * (1-mu[l])
  indH <- 1-indH; diag(indH) <- 0
  theta <- theta + indH * mu[l]

  beta.old <- log(matrix(1, N, N))

  g <- -log(rowSums(Pi*rho[l-1]*theta))

  while(l>t) { cat("B ", l, "\n")
    l <- l-1

    beta.new <- log(1.0 + (1-rho[l])*theta*exp(beta.old+g)) - g
    diag(beta.new) <- -g

    indH <- outer(H[,l], H[,l], "==") * 1; diag(indH) <- 0
    theta <- indH * (1-mu[l])
    indH <- 1-indH; diag(indH) <- 0
    theta <- theta + indH * mu[l]

    if(l>1) {
      g <- -( log(rowSums(Pi*rho[l-1]*theta*exp(beta.new+g))) - g)
    }

    beta.old <- beta.new
  }

  return(beta.old)
}
# Gold master forward algorithm, no sanity checks, matrix-wise, log-space
Exact_Forward_GMmatLS_R <- function(t, L, N, H, Pi, mu, rho) {
  # Setup Pi & rho
  diag(Pi) <- 0

  l <- 1

  indH <- outer(H[,l], H[,l], "==") * 1; diag(indH) <- 0
  theta <- indH * (1-mu[l])
  indH <- 1-indH; diag(indH) <- 0
  theta <- theta + indH * mu[l]

  alpha.old <- log(Pi*theta)

  f <- -log(rowSums(exp(alpha.old)*rho[l]))

  while(l<t) { cat("F ", l, "\n")
    l <- l+1

    indH <- outer(H[,l], H[,l], "==") * 1; diag(indH) <- 0
    theta <- indH * (1-mu[l])
    indH <- 1-indH; diag(indH) <- 0
    theta <- theta + indH * mu[l]

    alpha.new <- log(theta*Pi + theta*(1-rho[l-1])*exp(alpha.old+f)) - f

    f <- -( log(rowSums(exp(alpha.new+f)*rho[l])) - f )

    alpha.old <- alpha.new
  }

  return(list(alpha=alpha.old, alphasum=f))
}
