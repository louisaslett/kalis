#' Fast Calculation of Matrix Trace and Hilbert Schmidt Norm
#'
#' Provides multithreaded calculation of trace and Hilbert Schmidt Norm of a matrix PMP (where P is a projection matrix) without explicitly forming PMP.
#'
#' P here is assumed to have the form I-QQ' for some matrix Q of orthogonal columns
#'
#' @param M
#'        a real symmetric R matrix
#' @param tX
#'        `t((Q %*% (J%*%Q)) - (M %*% Q))`
#' @param tQ
#'        `t(Q)`
#' @param J
#'        `crossprod(Q, M)`
#' @param from_recipient
#'        haplotype index at which to start trace calculation -- useful for distributed computation (experimental feature, more documentation to come TODO)
#' @param nthreads
#'        the number of CPU cores to use.
#'        By default uses the `parallel` package to detect the number of physical cores.
#'
#' @return
#' a list containing three elements, the first is the trace `tr(PMP)`, the second is the *squared* Hilbert Schmidt Norm of PMP `tr((PMP)'PMP)`, the third is the diag of `PMP`.
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
