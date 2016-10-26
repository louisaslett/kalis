ComputeTable <- function(l, Pi, mu, rho, Ne = 20, method = "GM_LS_R") {
  # Transform rho
  rho <- as.matrix(rho)
  rho <- c(c(1-exp(-( (rho[2:L,1]-rho[1:(L-1),1]) * Ne * rho[1:(L-1),2] ))), 1)

  res <- switch(method,
    "GM_LS_R" = Exact_ComputeTable_GMmatLS_R(l, Pi, mu, rho),
    "C_naive" = Exact_ComputeTable_naive_C(l, Pi, mu, rho)
  )
  -(res$fwd+res$rev+res$alphasum)
}
