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
