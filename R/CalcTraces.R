#' Title TODO
#'
#' Short Description TODO
#'
#' Long Description TODO
#'
#' @param M
#'        TODO
#' @param tX
#'        TODO
#' @param tQ
#'        TODO
#' @param J
#'        TODO
#' @param from_recipient
#'        TODO
#' @param nthreads
#'        the number of CPU cores to use.
#'        By default uses the `parallel` package to detect the number of physical cores.
#'
#' @return
#' TODO
#'
#' @examples
#' # TODO
#'
#' @export
CalcTraces <- function(M, tX, tQ, J,
                       from_recipient = 1L,
                       nthreads = min(parallel::detectCores(logical = FALSE), ncol(M))) {
  .Call(CCall_CalcTraces, M, tX, tQ, J, from_recipient, nthreads)
}
