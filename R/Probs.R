#' Posterior marginal probabilities
#'
#' Calculate the posterior marginal probabilities at a given variant using forward and backward tables propagated to that position.
#'
#' The forward and backward tables must be at the same variant in order for them to be combined to yield the posterior marginal probabilities at variant \eqn{l}.
#' The \eqn{(j,i)}-th element of of the returned matrix is the probability that \eqn{j} is copied by \eqn{i} at the current variant, \eqn{l}, of the two tables, given the haplotypes observed (over the whole sequence).
#'
#' Note that each column represents an independent HMM.
#'
#' By convention, every diagonal element is zero.
#'
#' **Notes on `beta.theta.opts`**
#'
#' In order to obtain posterior marginal probability matrices between variants `fwd$l` and `bck$l`, then `bck` must be in "beta-theta space", see [Backward()] for details.
#' This allows the forward and backward tables to be transitioning both tables to some genomic position between `fwd$l` and `bck$l`.
#' The precise recombination distance by which each table is propagated can be controlled by passing optional arguments in a list via `beta.theta.opts`.
#' The recombination distances used can be specified in one of two ways.
#'
#' 1. Manually.
#'   In this case, `beta.theta.opts` is a list containing two named elements:
#'     - `"rho.fwd"` \eqn{\in (0,1)} specifies the transition probability `rho` for propagating the forward table.
#'     - `"rho.bck"` \eqn{\in (0,1)} specifies the transition probability `rho` for propagating the backward table.
#'
#' 2. Implicitly.
#'   In this case, `beta.theta.opts` is a list containing two named elements:
#'   - `"pars"`: a `kalisParameters` object that implicitly defines the recombination distance \eqn{\rho^\star} between `fwd$l` and `bck$l`
#'   - `"bias"` \eqn{\in (0,1)}.
#'     The forward table is propagated a distance of `(bias)`\eqn{\rho^\star} and the backward table is propagated a distance of `(1-bias)`\eqn{\rho^\star}.
#'
#' **Performance notes**
#'
#' If calculating many posterior probability matrices in succession, providing a pre-existing matrix `M` that can be updated in-place can dramatically increase speed by eliminating the time needed for memory allocation.
#' Be warned, since the matrix is updated in-place, if any other variables point to the same memory address, they will also be simultaneously overwritten.
#' For example, writing
#' ```
#' M <- matrix(0, nrow(fwd$alpha), ncol(fwd$alpha))
#' P <- M
#' PostProbs(fwd, bck, M = M)
#' ```
#' will update `M` and `P` simultaneously.
#'
#' When provided, `M` must have dimensions matching that of `fwd$alpha`.
#' Typically, that is simply \eqn{N \times N}{N x N} for \eqn{N} haplotypes.
#' However, if kalis is being run in a distributed manner, `M` will be a \eqn{N \times R}{N x R} matrix where \eqn{R} is the number of recipient haplotypes on the current machine.
#'
#' @param fwd a forward table as returned by [MakeForwardTable()] and propagated to a target variant by [Forward()].
#'   Must be at the same variant as `bck` (unless `bck` is in "beta-theta space" in which case if must be downstream ... see [Backward()] for details).
#' @param bck a backward table as returned by [MakeBackwardTable()] and propagated to a target variant by [Backward()].
#'   Must be at the same variant as `fwd` (unless `bck` is in "beta-theta space" in which case if must be downstream ... see [Backward()] for details).
#' @param unif.on.underflow a logical; if `TRUE`, then if all probabilities in a column underflow, then they will be set to \eqn{\frac{1}{N-1}}{1/(N-1)} instead of 0
#' @param M a pre-existing matrix into which to write the probabilities, can yield substantial speed up but requires special attention (see Details)
#' @param beta.theta.opts a list; see Details.
#' @param nthreads the number of CPU cores to use.
#'   By default uses the `parallel` package to detect the number of physical cores.
#'
#' @return
#'   Matrix of posterior marginal probabilities.
#'   The \eqn{(j,i)}-th element of of the returned matrix is the probability that \eqn{j} is copied by \eqn{i} at the current variant, \eqn{l}, of the two tables, given the haplotypes observed (over the whole sequence).
#'
#' @seealso
#'   [DistMat()] to generate calculate \eqn{d_{ji}}{d_(j,i)} distances directly;
#'   [Forward()] to propagate a Forward table to a new variant;
#'   [Backward()] to propagate a Backward table to a new variant
#'
#' @examples
#' # To get the posterior probabilities at, say, variants 100 of the toy data
#' # built into kalis
#' data("SmallHaps")
#' data("SmallMap")
#'
#' CacheHaplotypes(SmallHaps)
#'
#' rho <- CalcRho(diff(SmallMap))
#' pars <- Parameters(rho)
#'
#' fwd <- MakeForwardTable(pars)
#' bck <- MakeBackwardTable(pars)
#'
#' Forward(fwd, pars, 100)
#' Backward(bck, pars, 100)
#'
#' p <- PostProbs(fwd, bck)
#' d <- DistMat(fwd, bck)
#'
#' @export
PostProbs <- function(fwd, bck, unif.on.underflow = FALSE, M = NULL, beta.theta.opts = NULL,
                      nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient-fwd$from_recipient+1)){

  if(identical(nthreads, "R")) {
    if(!is.null(M)){stop("M cannot be NULL when requesting the gold master R version with R nthreads")}
    warning("Warning: using gold master R implementation.")
    return(invisible(PostProbs.GM(fwd, bck, unif.on.underflow, beta.theta.opts)))
  }

  rho.list <- input_checks_for_probs_and_dist_mat(fwd,bck,beta.theta.opts)

  # Make M if needed
  if(is.null(M)){ M <- matrix(0,nrow=nrow(fwd$alpha),ncol=ncol(fwd$alpha)) }

  # All clear to calculate distance matrices
  if(bck$beta.theta){
    vector.biproduct <- .Call(CCall_MatAndMulBtwVar, M, fwd, bck, rep(1, nrow(M)), FALSE, FALSE, TRUE, unif.on.underflow, rho.list$rho.fwd, rho.list$rho.bck, nthreads)
  } else {
    vector.biproduct <- .Call(CCall_MatAndMul, M, fwd, bck, rep(1, nrow(M)), FALSE, FALSE, TRUE, unif.on.underflow, nthreads)
  }

  invisible(M)
}

#' Distance matrix
#'
#' Utility for calculating distance matrices at, in between, or excluding variants.
#'
#' This computes a local probability or distance matrix based on the forward and backward tables at a certain variant.
#' The default usage is provide forward and backward tables at the same variant \eqn{l} so that the \eqn{(j,i)}-th element of of the returned matrix is the inferred distance \eqn{d_{ji}}{d_(j,i)} between haplotypes \eqn{j} and \eqn{i} at the current variant, \eqn{l}, of the two tables given the haplotypes observed (over the whole sequence).
#'
#' In particular,
#'
#' \deqn{d_{ji} = -log(p_{ji})}{d_(j,i) = - log(p_(j,i)) }
#'
#' where \eqn{p_{ji}}{p_(j,i)} is the posterior marginal probability that \eqn{j} is coped by \eqn{i} at the current variant of the two tables, \eqn{l}, given the haplotypes observed (over the whole sequence).
#'
#' By convention, \eqn{d_{ii} = 0}{d_(i,i) = 0} for all \eqn{i}.
#'
#' The above is returned when the `type` argument is `"raw"` (the default).
#' However, for convenience, a user may change the `type` argument to `"std"` in order for the distances to be mean and variance normalized before returning.
#' Changing `type` to `"minus.min"` will subtract the min distance in each column from that column before returning (this is the default in RELATE, see bottom of the RELATE paper Supplement Page 7 (Speidel et al, 2019)).
#'
#' This function also allows users to calculate distance matrices in between variants and also to calculate matrices that exclude a set of consecutive variants by passing a backward table in "beta-theta space."
#' If in "beta-theta space", `bck$l` may be greater than but not equal to `fwd$l`.
#' `beta.theta.opts` provides is required in this case to set how much of a recombination distance to propagate each matrix before combining them into distances.  See Details below.
#'
#' **Notes on `beta.theta.opts`**
#'
#' In order to obtain distance matrices between variants `fwd$l` and `bck$l`, then `bck` must be in "beta-theta space", see [Backward()] for details.
#' This allows the forward and backward tables to be transitioning both tables to some genomic position between `fwd$l` and `bck$l`.
#' The precise recombination distance by which each table is propagated can be controlled by passing optional arguments in a list via `beta.theta.opts`.
#' The recombination distances used can be specified in one of two ways.
#'
#' 1. Manually.
#'   In this case, `beta.theta.opts` is a list containing two named elements:
#'     - `"rho.fwd"` \eqn{\in (0,1)} specifies the transition probability `rho` for propagating the forward table.
#'     - `"rho.bck"` \eqn{\in (0,1)} specifies the transition probability `rho` for propagating the backward table.
#'
#' 2. Implicitly.
#'   In this case, `beta.theta.opts` is a list containing two named elements:
#'   - `"pars"`: a `kalisParameters` object that implicitly defines the recombination distance \eqn{\rho^\star} between `fwd$l` and `bck$l`
#'   - `"bias"` \eqn{\in (0,1)}.
#'     The forward table is propagated a distance of `(bias)`\eqn{\rho^\star} and the backward table is propagated a distance of `(1-bias)`\eqn{\rho^\star}.
#'
#' **Performance notes**
#'
#' If calculating many posterior probability matrices in succession, providing a pre-existing matrix `M` that can be updated in-place can dramatically increase speed by eliminating the time needed for memory allocation.
#' Be warned, since the matrix is updated in-place, if any other variables point to the same memory address, they will also be simultaneously overwritten.
#' For example, writing
#' ```
#' M <- matrix(0, nrow(fwd$alpha), ncol(fwd$alpha))
#' P <- M
#' PostProbs(fwd, bck, M = M)
#' ```
#' will update `M` and `P` simultaneously.
#'
#' When provided, `M` must have dimensions matching that of `fwd$alpha`.
#' Typically, that is simply \eqn{N \times N}{N x N} for \eqn{N} haplotypes.
#' However, if kalis is being run in a distributed manner, `M` will be a \eqn{N \times R}{N x R} matrix where \eqn{R} is the number of recipient haplotypes on the current machine.
#'
#' @param fwd a forward table as returned by [MakeForwardTable()] and propagated to a target variant by [Forward()].
#'   Must be at the same variant as `bck` (unless `bck` is in "beta-theta space" in which case if must be downstream ... see [Backward()] for details).
#' @param bck a backward table as returned by [MakeBackwardTable()] and propagated to a target variant by [Backward()].
#'   Must be at the same variant as `fwd` (unless `bck` is in "beta-theta space" in which case if must be downstream ... see [Backward()] for details).
#' @param type a string; must be one of `"raw"`, `"std"` or `"minus.min"`.
#'   See Details.
#' @param M a pre-existing matrix into which to write the distances.
#'   This can yield substantial speed up but requires special attention, see Details.
#' @param beta.theta.opts a list; see Details.
#' @param nthreads the number of CPU cores to use. By default no parallelism is used.
#'   By default uses the `parallel` package to detect the number of physical cores.
#'
#' @return A matrix of distances.
#'   The \eqn{(j,i)}-th element of of the returned matrix is the inferred
#'   distance \eqn{d_{ji}}{d_(j,i)} to haplotype \eqn{j} from haplotype \eqn{i} at the current variant.
#'   Each column encodes the output of an independent HMM: in column \eqn{i}, haplotype \eqn{i} is taken as the observed recipient haplotype and painted as a mosaic of the other \eqn{N-1} haplotypes.
#'   Hence, the distances are asymmetric.
#'
#'   If you wish to plot this matrix or perform clustering, you may want to symmetrize the matrix first.
#'
#' @references
#'   Speidel, L., Forest, M., Shi, S., & Myers, S. (2019). A method for genome-wide genealogy estimation for thousands of samples. *Nature Genetics*, **51**(1321–1329).
#'
#' @seealso
#'   [PostProbs()] to calculate the posterior marginal probabilities \eqn{p_{ji}}{p_(j,i)};
#'   [Forward()] to propagate a Forward table to a new variant;
#'   [Backward()] to propagate a Backward table to a new variant.
#'
#' @examples
#' # To get the posterior probabilities at, say, variants 100 of the toy data
#' # built into kalis
#' data("SmallHaps")
#' data("SmallMap")
#'
#' CacheHaplotypes(SmallHaps)
#'
#' rho <- CalcRho(diff(SmallMap))
#' pars <- Parameters(rho)
#'
#' fwd <- MakeForwardTable(pars)
#' bck <- MakeBackwardTable(pars)
#'
#' Forward(fwd, pars, 100)
#' Backward(bck, pars, 100)
#'
#' p <- PostProbs(fwd, bck)
#' d <- DistMat(fwd, bck)
#'
#' \donttest{
#' plot(d)
#' }
#'
#' @export DistMat
DistMat <- function(fwd, bck, type = "raw", M = NULL, beta.theta.opts = NULL,
                    nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient-fwd$from_recipient+1)){

  if(identical(nthreads, "R")) {
    if(!is.null(M)){stop("M must be NULL when requesting the gold master R version with R nthreads")}
    warning("Warning: using gold master R implementation.")
    return(invisible(DistMat.GM(fwd, bck, type, beta.theta.opts)))
  }

  rho.list <- input_checks_for_probs_and_dist_mat(fwd,bck,beta.theta.opts)

  # Make M if needed
  if(is.null(M)){ M <- matrix(0,nrow=nrow(fwd$alpha),ncol=ncol(fwd$alpha)) }

  if(type == "raw") {
    stdz <- FALSE
    minus.min <- FALSE
  } else if(type == "std") {
    stdz <- TRUE
    minus.min <- FALSE
  } else if(type == "minus.min") {
    stdz <- FALSE
    minus.min <- TRUE
  } else {
    stop(glue("argument type '{type}' is not a recognised option"))
  }

  # All clear to calculate distance matrices
  if(bck$beta.theta){
    vector.biproduct <- .Call(CCall_MatAndMulBtwVar, M, fwd, bck, rep(1, nrow(M)), stdz, minus.min, FALSE, FALSE, rho.list$rho.fwd, rho.list$rho.bck, nthreads)
  } else {
    vector.biproduct <- .Call(CCall_MatAndMul, M, fwd, bck, rep(1, nrow(M)), stdz, minus.min, FALSE, FALSE, nthreads)
  }

  invisible(M)
}



input_checks_for_probs_and_dist_mat <-  function(fwd,bck,beta.theta.opts = NULL){

  # RUN GENERAL CHECKS
  if(fwd$l == 2147483647L){stop("forward table has not been initialized but not propagated to a variant in {1,...,L}.")}
  if(bck$l == 2147483647L){stop("backward table has not been initialized but not propagated to a variant in {1,...,L}.")}

  if(fwd$pars.sha256 != bck$pars.sha256) {
    stop("parameters used to calculate the forward table and backward table do not match.")
  }
  if(fwd$from_recipient != bck$from_recipient || fwd$to_recipient != bck$to_recipient) {
    stop("forward and backward tables cover different from/to recipient windows.")
  }

  if(fwd$l > bck$l){stop("fwd$l > bck$l.  The forward table cannot be past the backward table.")}


  # RUN BTW LOCI and AT LOCI specific checks
  if(bck$beta.theta){

    if(fwd$l == bck$l){stop("A forward table cannot be combined with a backward table at the same variant if the backward table is in beta.theta space.")}

    if(!((!is.null(beta.theta.opts$rho.fwd) & !is.null(beta.theta.opts$rho.bck)) | (!is.null(beta.theta.opts$pars) & !is.null(beta.theta.opts$bias)))){
      stop("beta.theta.opts must be a named list containing either pars and bias OR rho.fwd and rho.bck.")
    }

    if(!is.null(beta.theta.opts$rho.fwd) & !is.null(beta.theta.opts$rho.bck)){
      # Use custom morgan distances

      if(!is.numeric(beta.theta.opts$rho.fwd) || beta.theta.opts$rho.fwd<=0 ){stop("rho.fwd must be numeric and strictly positive.  To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}
      if(!is.numeric(beta.theta.opts$rho.bck) || beta.theta.opts$rho.bck<=0 ){stop("rho.bck must be numeric and strictly positive.  To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}

      rho.fwd <- beta.theta.opts$rho.fwd
      rho.bck <- beta.theta.opts$rho.bck

    }else{
      # Use bias and pars

      if(!inherits(beta.theta.opts$pars,"kalisParameters")){stop("beta.theta.opts$pars must be kalisParameters object.")}

      if(!is.numeric(beta.theta.opts$bias) || beta.theta.opts$bias<0 || beta.theta.opts$bias>1 ){stop("bias must be numeric and within [0,1]. To obtain a distance matrix AT a particular variant, advance bck to that variant in beta space.")}

      total.rho <- sum(beta.theta.opts$pars$pars$rho[fwd$l:(bck$l - 1)])

      rho.fwd <- total.rho * beta.theta.opts$bias
      rho.bck <- total.rho * (1-beta.theta.opts$bias)

    }

    return(invisible(list("rho.fwd" = rho.fwd, "rho.bck" = rho.bck)))

  }else{

    if(bck$l != fwd$l){stop("variant position of the forward table and backward table do not match.")}
    return(invisible(NULL))
  }
}



#' Plotting function for a kalisDistanceMatrix object
#'
#' Clusters the given distance matrix and generates a heatmap to display it.
#'
#' @param d a kalisDistanceMatrix
#'
#' @return There is nothing returned.
#'
#' @export
plot.kalisDistanceMatrix <- function(x, cluster.method = "average", ...) {

  perm <- fastcluster::hclust(stats::as.dist(x),method=cluster.method)$order
  print(lattice::levelplot(x[perm,][,rev(perm)],
                           useRaster = TRUE,
                           col.regions = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,name = "BuPu"))(100),
                           yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxt = "n", ...))
}
