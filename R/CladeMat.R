#' Fast Clade Matrix Construction
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
