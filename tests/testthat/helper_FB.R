start.new.example <- function(N, L, mut.rate = 0.3, Ne = 1, gamma = 1) {
  suppressWarnings(CacheHaplotypes(matrix(sample(0:1, N*L, replace = TRUE), nrow = L, ncol = N)))

  # This is a not too bad fit to the empirical distribution of Morgan distances
  # for some real recombination maps, just to get something approximately
  # realistic for testing (NB this should really change with the length of
  # sequence)
  Morgan.dist <- exp(rlogis(L-1, -8.8384, 1.0920))
  while(max(Morgan.dist) > 0.008)
    Morgan.dist[Morgan.dist > 0.008] <- exp(rlogis(sum(Morgan.dist > 0.008), -8.8384, 1.0920))
  rho <- CalcRho(cM = 100.0*Morgan.dist, s = Ne, gamma = gamma)

  mu <- rnorm(L, 1e-8, 1e-8)
  mu[1] <- 1e-8

  Pi <- matrix(abs(rnorm(N*N, 1/(N-1), 1/(10*N-1))), nrow = N, ncol = N)
  attr.Pi <- attributes(Pi)
  diag(Pi) <- 0
  Pi <- scale(Pi, center = FALSE, scale = colSums(Pi))
  attributes(Pi) <- attr.Pi

  list(rho = rho, mu = mu, Pi = Pi)
}

make.pars <- expression({
  pars <- Parameters(p$rho, mu = p$mu, Pi = p$Pi, use.speidel = use.speidel)
  pars.scmu <- Parameters(p$rho, mu = p$mu[1], Pi = p$Pi, use.speidel = use.speidel)
  pars.scPi <- Parameters(p$rho, mu = p$mu, Pi = NULL, use.speidel = use.speidel)
  pars.scmuPi <- Parameters(p$rho, mu = p$mu[1], Pi = NULL, use.speidel = use.speidel)
})

err <- function(aprx, trth) {
  if(!identical(dim(aprx), dim(trth))) {
    stop("Dimensions don't match")
  }
  abs.err <- abs(c(trth-aprx))
  rel.err <- abs.err[c(trth)!=0]/c(trth[c(trth)!=0])
  list(abs.err = max(abs.err), rel.err = max(rel.err), zero.mismatches = sum(aprx[c(trth)==0]!=0))
}
