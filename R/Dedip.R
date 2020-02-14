#' Dedip in Parallel
#' @export
DedipPar <- function(M, x = rep(1/nrow(M[[1]]),nrow(M[[1]])), fwd, bck, pars, fwd_rho = NULL, bck_rho = NULL, from_recipient = 1, nthreads = 1){
  # fwd.rho and bck.rho allow you to calculate rho corresponding to whatever
  # cM distance you want to go so that you can easily go to an exact cM coordinate if you wish.

  if(length(M)!=4){stop("M must be a list of four matrices to hold the dediploided results.")}

  if(fwd$l == bck$l){
    if(bck$beta.theta){stop("Cannot calculate haplotype distances when fwd$l == bck$l if bck is in beta.theta space : that would double-count the variant at l.")}

    DedipParAtVarScalarPi(M, fwd, bck, x, from_recipient, nthreads)

  }else{

    if(fwd$l > bck$l){stop("fwd$l is > bck$l : the forward table cannot be past the backward table.")}
    if(class(pars$pars$Pi)=="matrix"){stop("The methods for combining fwd and bck when fwd$l < bck$l have not been implemented for the case where Pi is a general matrix.")}
    if(!bck$beta.theta){stop("The backward table must be in beta.theta space for cases where fwd$l < bck$l.")}

    if(is.null(fwd_rho) & is.null(bck_rho)){
      # Calculate the probability of at least one transition between fwd$l and bck$l and divide it by 2 to get fwd_rho == bck_rho.
      fwd_rho <- bck_rho <- -expm1(sum(log1p(-pars$pars$rho[fwd$l:(bck$l-1)]))) / 2
    } else {
      if(is.null(fwd_rho) | is.null(bck_rho)){stop("If fwd_rho or bck_rho is specified than the other must be provided also.")}
      if(!is.numeric(fwd_rho) | !is.numeric(bck_rho)){stop("Both fwd_rho and bck_rho must be numeric.")}
      if(!is.finite(fwd_rho) | !is.finite(bck_rho)){stop("Both fwd_rho and bck_rho must be finite.")}
    }

    DedipParBtwVarScalarPi(M, fwd, bck, x, fwd_rho, bck_rho, from_recipient, nthreads)

  }
}






Dedip <- function(fwd, bck, method="min") {

  s <- colSums(fwd$alpha*bck$beta)

  return(switch(method,
                min    = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=Dedip_min(fwd$alpha, bck$beta, s), uplo="L"),
                min2nd = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=Dedip_2nd_min(fwd$alpha, bck$beta, s), uplo="L"),
                dom    = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=Dedip_dom(fwd$alpha, bck$beta, s), uplo="L"),
                max    = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=Dedip_max(fwd$alpha, bck$beta, s), uplo="L"),
                add    = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=Dedip_add(fwd$alpha, bck$beta, s), uplo="L"),
                mean   = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=Dedip_mean(fwd$alpha, bck$beta, s), uplo="L"),
                all    = {
                  res <- Dedip_all(fwd$alpha, bck$beta, s)
                  list(min    = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=res[[1]], uplo="L"),
                       min2nd = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=res[[2]], uplo="L"),
                       dom    = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=res[[3]], uplo="L"),
                       max    = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=res[[4]], uplo="L"),
                       add    = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=res[[5]], uplo="L"),
                       mean   = methods::new("dspMatrix", Dim=as.integer(c(nrow(fwd$alpha)/2,nrow(fwd$alpha)/2)), x=res[[6]], uplo="L"))
                },
                stop("Unknown method")))
}

Dedip2 <- function(M, method="min") {
# input to Dedip2 is a matrix ( -log10(alpha * beta) - mu )/ sigma that's already been column standardized
  return(switch(method,
                min    = methods::new("dspMatrix", Dim=as.integer(c(nrow(M)/2,nrow(M)/2)), x=Dedip2_min(M), uplo="L"),
                min2nd = methods::new("dspMatrix", Dim=as.integer(c(nrow(M)/2,nrow(M)/2)), x=Dedip2_2nd_min(M), uplo="L"),
                dom    = methods::new("dspMatrix", Dim=as.integer(c(nrow(M)/2,nrow(M)/2)), x=Dedip2_dom(M), uplo="L"),
                all    = {
                  res <- Dedip2_all(M)
                  list(min    = methods::new("dspMatrix", Dim=as.integer(c(nrow(M)/2,nrow(M)/2)), x=res[[1]], uplo="L"),
                       min2nd = methods::new("dspMatrix", Dim=as.integer(c(nrow(M)/2,nrow(M)/2)), x=res[[2]], uplo="L"),
                       dom    = methods::new("dspMatrix", Dim=as.integer(c(nrow(M)/2,nrow(M)/2)), x=res[[3]], uplo="L"))
                },
                stop("Unknown method")))
}