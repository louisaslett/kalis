#' @export
CalcTraces <- function(M, tX, tQ, J,
                       from_recipient = 1L,
                       nthreads = min(parallel::detectCores(logical = FALSE), ncol(M))) {
  .Call(CCall_CalcTraces, M, tX, tQ, J, from_recipient, nthreads)
}
