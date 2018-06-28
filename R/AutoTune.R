#' Auto-tune the parameters for HMM
#'
#' Provides auto-tuning on a given chromosome for the \code{Ne}, \code{mu} and
#' \code{gamma} parameters of the HMM.
#'
#' Detailed description
#'
#' Also reference the vignette
#'
#' @param t the number of loci to auto-tune against (which will be evently spaced
#'   according to the inverse cumulative recombination map); or a vector of specific
#'   loci to use for auto-tuning.
#' @param cache a forward table cache as returned by \code{\link{CreateForwardTableCache}}.
#'   Note this cache will be overwritten.
#' @param morgan.dist a vector of recombination distances between loci, in Morgans.
#'   Note element i of this vector should be the distance between loci i and i+1
#'   (not i and i-1), and thus length one less than the sequence length.
#' @param Pi leaving the default of uniform copying probabilities is recommended for
#'   computational efficiency.  If desired, a full matrix of background copying
#'   probabilities can be provided, such that the (i,j)-th element is the background
#'   probability that i copies j.  Hence, (a) the diagonal must be zero; and (b)
#'   the columns of Pi must sum to 1.
#' @param nthreads the number of CPU cores on which to run.
#'
#' @return The optimal parameters for the HMM as a list.
#'
#' @seealso \code{\link{Forward}}, \code{\link{ForwardUsingTableCache}}, and
#'   \code{\link{Backward}} for functions which consume the parameters which
#'   can be auto-tuned by this function.
#'
#' @examples
#'
AutoTune <- function(t, cache, morgan.dist, Pi = 1/(nrow(fwd$alpha)-1), nthreads = 1) {
  if(!is.vector(t)) {
    stop("t must be either a vector or a scalar.")
  }
  if(!is.numeric(t)) {
    stop("t must be numeric.")
  }
  if(length(t) < 1 || length(t) > L) {
    stop("t is the wrong length for this problem.")
  }


  if(length(t) == 1) {
    t <-
  }

}

InvRecombMap <- function(morgan.dist,num.loci) {
  #morgan.dist<-c(0.1,0.4,0.2,0.23,0.43,0.12,0.7,0.10)
  morgan.cdf<-cumsum(morgan.dist)
  
  # Could just set it so that if 
  
  #num.loci<-6
  
  # If an two SNPs are further apart than the gap between target loci, then we don't want to depend 
  # on the estimation of that recombination hotspot, so we simply ignore that hotspot to be convervative
  
  delta<-(morgan.cdf[length(morgan.cdf)]-morgan.cdf[1])/(num.loci+1)
  # loci.to.exclude<-which(morgan.dist>(delta-1e-12))
  # morgan.dist.thinned<-morgan.dist[-loci.to.exclude]
  # locs<-2:(length(morgan.dist)+1)[-loci.to.exclude]
  # morgan.cdf.thinned<-cumsum(morgan.dist.thinned)
  # 
  target.grid<-seq(0,morgan.cdf[length(morgan.cdf)],length.out = num.loci+2)[-1]
  target.grid<-target.grid[-length(target.grid)]

  yy<-approxfun(x=morgan.cdf,y=2:(length(morgan.dist)+1),method="constant")(target.grid)
  #plot(2:(length(morgan.dist)+1),morgan.cdf,type="S")
  #for(i in 1:length(target.grid)){abline(h=target.grid[i])}
  
  #for(i in 1:length(yy)){abline(v=yy[i])}
  return(yy)
  # ties 
  # integer
  # robust to all possible inputs
  }
