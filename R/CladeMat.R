#' Fast clade matrix construction
#'
#' Constructs a clade matrix using forward and backward tables.
#' The clade matrix captures genetic relatedness information in the distances from the Li & Stephens model that are not captured in the called clades.
#'
#' `CladeMat()` uses the forward and backward tables to construct the corresponding clade matrix which can then be tested, for example using a standard quadratic form score statistic.
#'
#' @references
#' Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I. (2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv 2024.09.30.615852. Available at: \doi{10.1101/2024.09.30.615852}.
#'
#' @param fwd
#'        a `kalisForwardTable` object, as returned by [MakeForwardTable()] and propagated to a target variant by [Forward()].
#'        This table must be at the same variant location as argument `bck`.
#' @param bck
#'        a `kalisBackwardTable` object, as returned by [MakeBackwardTable()] and propagated to a target variant by [Backward()].
#'        This table must be at the same variant location as argument `fwd`.
#' @param M
#'        a matrix with half the number of rows and columns as the corresponding forward/backward tables.
#'        This matrix is overwritten in place with the clade matrix result for performance reasons.
#' @param unit.dist
#'        the change in distance that is expected to correspond to a single mutation (typically \eqn{-\log(\mu)}) for the LS model)
#' @param thresh
#'        a regularization parameter: differences of distances must exceed this threshold (in `unit.dist` units) in order to cause the introduction of a probabilistic clade.
#'        Defaults to `0.2`.
#' @param max1var
#'        a logical regularization parameter.
#'        When `TRUE`, differences in distances exceeding 1 `unit.dist` are set to 1 (so that any edge in the latent ancestral tree with multiple mutations on them are treated as if only one mutation was on it).
#' @param nthreads
#'        the number of CPU cores to use.
#'        By default uses the `parallel` package to detect the number of physical cores.
#'
#' @return
#' A list, the first element contains a list of tied nearest neighbours (one for each haplotype).
#' Other elements of the returned list are for internal use by [PruneCladeMat()] to allow for efficient removal of singletons and sprigs.
#'
#' @examples
#' # TODO
#'
#'
#' @export CladeMat
CladeMat <- function(fwd, bck, M, unit.dist, thresh = 0.2, max1var = FALSE,
                     nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient-fwd$from_recipient+1)){

  # input checks
  #########################
  input_checks_for_probs_and_dist_mat(fwd,bck)

  if(nrow(fwd$alpha)%%2 !=0 || ncol(fwd$alpha)%%2 !=0 || nrow(bck$beta)%%2 !=0 || ncol(bck$beta)%%2 !=0 ){
    stop("fwd and bck must both have an even number of recipient haplotypes and an even number of donor haplotypes")
  }

  if(!is.matrix(M) || !is.double(M) || nrow(M) != nrow(fwd$alpha)/2 || ncol(M) != ncol(fwd$alpha)/2){
    stop("M must be a matrix of doubles with nrow(fwd$alpha)/2 rows and ncol(fwd$alpha)/2 columns")}

  if(!is.atomic(unit.dist) || length(unit.dist)!=1L || !is.finite(unit.dist) || unit.dist <= 0){
    stop("unit.dist must be a number greater than 0")}

  if(is.integer(unit.dist)){
    unit.dist <- as.double(unit.dist)
  } else {
    if(!is.double(unit.dist)){stop("unit.dist must be a number greater than 0")}}

  if(!is.atomic(thresh) || length(thresh)!=1L || !is.finite(thresh) || thresh < 0 || thresh > 1){
    stop("thresh must be a number in [0,1]")}

  if(is.integer(thresh)){
    thresh <- as.double(thresh)
  } else {
    if(!is.double(thresh)){stop("thresh must be a number in [0,1]")}}

  if(!is.logical(max1var) || length(max1var) > 1){
    stop("max1var must be TRUE or FALSE")}

  nthreads <- as.integer(nthreads)
  if(!is.integer(nthreads) || length(nthreads)!=1L || !is.finite(nthreads) || nthreads <= 0){
    stop("nthreads must be a positive integer")}

  if(nthreads > ncol(fwd$alpha)/2){
    stop("nthreads cannot be greater than the number of recipient haplotypes divided by 2.")
  }

  invisible(.Call(CCall_CladeMat, fwd, bck, M, unit.dist, thresh, max1var, nthreads))
}
