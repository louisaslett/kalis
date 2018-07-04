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
AutoTune <- function(t, cache, morgan.dist, gamma = NULL, it = 30, nthreads = 1) {
  init_points <- 4

  if(!is.vector(t)) {
    stop("t must be either a vector or a scalar.")
  }
  if(!is.numeric(t)) {
    stop("t must be numeric.")
  }
  if(length(t) < 1 || length(t) >= L) {
    stop("t is the wrong length for this problem.")
  }
  if(length(t) == 1 && (t < 1 || t >= L)) {
    stop("Invalid number of loci chosen for auto-tuning.")
  }
  if(length(t) > 1 && !isTRUE(all.equal(t, sort(t)))) {
    stop("Loci specified in t are not sorted ... are these definitely loci for auto-tuning?")
  }
  if(length(t) > 1 && (any(t<1) || any(t>L))) {
    stop("Invalid loci numbers specified in t.")
  }
  if(class(cache) != "StatGenCheckpointTable") {
    stop("cache must be a StatGenCheckpointTable.")
  }
  if(!is.vector(morgan.dist)) {
    stop("morgan.dist must be a vector of recombination distances.")
  }
  if(!is.numeric(morgan.dist)) {
    stop("morgan.dist must be numeric vector type.")
  }
  if(length(morgan.dist) != L-1) {
    stop("morgan.dist is the wrong length for this problem.")
  }
  if(!is.null(gamma)) {
    if(!is.vector(gamma) || !is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
      stop("gamma must be a positive scalar.")
    }
  }
  # if(is.data.frame(Pi)) {
  #   stop("Pi must be a matrix or scalar, not a data frame.")
  # }
  # # if(is.matrix(Pi) && (nrow(Pi) != N || ncol(Pi) != N)) {
  # #   stop("Pi is of the wrong dimensions for this problem.")
  # if(is.matrix(Pi)) {
  #   stop("Only scalar Pi is supported for auto-tuning.")
  # }
  # if(!is.matrix(Pi) && !(is.vector(Pi) && is.numeric(Pi) && length(Pi) == 1 && Pi == 1/(nrow(fwd$alpha)-1))) {
  #   stop("Pi can only be set to a matrix, or omitted to have uniform copying probabilities of 1/(N-1) for a problem with N recipients.")
  # }
  Pi <- 1/(nrow(cache[[1]]$alpha)-1)

  # If not given specific loci, use the cumulative recombination map to find sensibly spaced ones
  if(length(t) == 1) {
    t <- InvRecombMap(morgan.dist, t)
  }

  # Check the cache is big enough -- we should add checkpointing in future?
  if(length(cache) < length(t)) {
    stop("cache is not large enough to store the forward table at all loci: checkpointing will be implemented in a future release, for now please increase table cache size or auto-tune using fewer loci.")
  }

  # Tuning combinations
  if(is.null(gamma)) {
    BO.fn <- function(Ne, gamma, mu) {
      print(c("Trying Ne =", Ne, ", gamma =", gamma, ", mu =", mu, "\n"))
      AutoTuneTarget(bck, cache, t, morgan.dist, Ne, gamma, mu, Pi, e, nthreads)
    }
    BO.bounds <- list(Ne = c(0.01, 100), gamma = c(0.1, 10), mu = c(1e-10, 1))
  } else {
    BO.fn <- function(Ne, mu) {
      print(c("Trying Ne =", Ne, ", mu =", mu, "\n"))
      AutoTuneTarget(bck, cache, t, morgan.dist, Ne, gamma, mu, Pi, e, nthreads)
    }
    BO.bounds <- list(Ne = c(0.01, 100), mu = c(1e-10, 1))
  }

  # We want to extract the negative log p-values for diagnostics
  e <- new.env(parent = emptyenv())
  e$neglog.p.val <- matrix(0, init_points+it, length(t))
  e$i <- 1

  # Do Bayesian optimisation search for best parameters
  bck <- MakeBackwardTable(cache[[1]]$from_recipient, cache[[1]]$to_recipient)
  bo <- BayesianOptimization(BO.fn, bounds = BO.bounds,
                             init_points = init_points, n_iter = it)
  rm(bck)
  list(BayesianOptimization = bo,
       neglog.p.val = e$neglog.p.val)
}

AutoTuneTarget <- function(bck, cache, t, morgan.dist, Ne, gamma, mu, Pi, e, nthreads) {
  neglog.p.val <- rep(0, length(t))

  rho <- c(1-exp(-Ne*morgan.dist^gamma), 1)
  rho <- ifelse(rho<1e-16, 1e-16, rho)

  ResetBackwardTable(bck)
  AutoTuneFillForwardCache(cache, t, morgan.dist, Ne, gamma, mu, Pi, nthreads)
  for(i in length(t):1) {
    Forward1step_scalarPi_scalarmu_cpp(cache[[i]], t[i], Pi, mu, rho, nthreads)
    Backward(bck, t[i], morgan.dist, Ne, gamma, mu, Pi, nthreads)

    # Get distance matrix
    M <- -log(t(cache[[i]]$alpha*bck$beta)) + cache[[i]]$alpha.f2 + bck$beta.g2
    diag(M) <- Inf
    min_dist_vec <- M[cbind(1:nrow(M), max.col(-M))]
    lZ <- log(rowSums(exp(min_dist_vec - M))) - min_dist_vec
    M <- M+lZ
    diag(M) <- 0
    M <- (M+t(M))/2

    # Compute p-value
    y <- QueryCache(start = t[i], length = 1)[,1]
    mu <- mean(y)
    sigma <- sqrt(mu*(1-mu))

    obs <- c(crossprod(y-mu, M%*%(y-mu)))
    neglog.p.val[i] <- -QFIntBounds(obs, as(M, "dsyMatrix"), rep(mu, nrow(M)), rep(sigma, nrow(M)), 25, log = TRUE)[1,4]
  }

  e$neglog.p.val[e$i,] <- neglog.p.val
  e$i <- e$i+1

  list(Score = sum(neglog.p.val), Pred = 0)
}

AutoTuneFillForwardCache <- function(cache, t, morgan.dist, Ne, gamma, mu, Pi, nthreads) {
  t.pre <- t-1

  # Use the final element of cache as the table to move forward
  N <- length(t)
  ResetForwardTable(cache[[N]])

  for(i in 1:(N-1)) {
    Forward(cache[[N]], t.pre[i], morgan.dist, Ne, gamma, mu, Pi, nthreads)
    CopyForwardTable(cache[[i]], cache[[N]])
  }
}

InvRecombMap <- function(morgan.dist, num.target.loci) {
  # USE BELOW COMMENTED CODE FOR EXAMPLE
  #morgan.dist.orig<-c(0.1,0.4,0.2,0.23,0.43,0.12,0.7,0.10, 0.12,0.14,0.02,0.23,1,0.22,1.2,0.9)
  #raw_map <- as.matrix(read.table(paste("Documents/lab/software/stat_gen/MalariaGEN_data/maps/genetic_map_chr2_combined_b37.txt"),header=T))
  #morgan.dist.orig<-diff(raw_map[,3])[1:1000]
  #morgan.dist.orig<-c(morgan.dist.orig[1:200],.5,morgan.dist.orig[201:600],.5,morgan.dist.orig[601:1000])
  #morgan.dist<-morgan.dist.orig
  #num.target.loci<-10

  if(length(morgan.dist) < num.target.loci) stop("Cannot return more target loci than exist in morgan.dist")

  L <- length(morgan.dist)+1
  delta <- (sum(morgan.dist)-morgan.dist[1])/(num.target.loci+1)-sqrt(.Machine$double.eps)

  # This function finds all points in the recombination map where the map jumps more than delta and
  # then tempers the peaks of the recomb map just enough to ensure that all of the peaks are less than delta
  # while redistributing the mass so that the recombination map CDF adds up to the same number.
  # This allows us to maintain the local recombination architecture while ensuring that training points are all
  # distinct integers that do not neglect sampling more densely from regions with high recombination hotspots

  if(any(morgan.dist>delta)) {
    tryCatch({
      ideal.x <- uniroot(function(x, morgan.dist, L, delta) {
        temp <- (morgan.dist-x)
        eps <- sum(temp[which(temp>0)])/(L-1)
        return(delta-(x+eps))
      }, morgan.dist = morgan.dist, L = L, delta = delta, upper = delta, lower = 0)$root
    }, error = function(e) stop("Cannot uniformly space that many target loci on this morgan.dist object.  Try reducing num.target.loci."))
    temp <- morgan.dist-ideal.x
    eps <- sum(temp[which(temp>0)])/(L-1)
    morgan.dist <- ifelse(morgan.dist<=ideal.x, morgan.dist, ideal.x) + eps
  }

  morgan.cdf <- cumsum(morgan.dist)
  target.grid <- seq(0,morgan.cdf[L-1], length.out = num.target.loci+2)[-1]
  target.grid <- target.grid[-length(target.grid)]

  yy <- approxfun(x = morgan.cdf, y = 2:L, method = "constant")(target.grid)

  if(anyDuplicated(yy)!=0) stop("Duplicated target loci returned when inverting recomb map")

  #plot(2:L,cumsum(morgan.dist.orig),type="S")
  #lines(2:L,morgan.cdf,type="S",col="blue")
  # for(i in 1:length(target.grid)){abline(h=target.grid[i])}
  # for(i in 1:length(yy)){abline(v=yy[i])}

  yy
}
