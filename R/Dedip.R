Dedip <- function(fwd, bck, method="min") {
  return(switch(method,
                min    = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=Dedip_min(fwd, bck), uplo="L"),
                min2nd = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=Dedip_2nd_min(fwd, bck), uplo="L"),
                dom    = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=Dedip_dom(fwd, bck), uplo="L"),
                max    = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=Dedip_max(fwd, bck), uplo="L"),
                add    = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=Dedip_add(fwd, bck), uplo="L"),
                mean   = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=Dedip_mean(fwd, bck), uplo="L"),
                all    = {
                  res <- Dedip_all(fwd, bck)
                  list(min    = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=res[[1]], uplo="L"),
                       min2nd = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=res[[2]], uplo="L"),
                       dom    = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=res[[3]], uplo="L"),
                       max    = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=res[[4]], uplo="L"),
                       add    = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=res[[5]], uplo="L"))
                       mean   = new("dspMatrix", Dim=as.integer(c(nrow(fwd)/2,nrow(fwd)/2)), x=res[[6]], uplo="L"))
                },
                stop("Unknown method")))
}