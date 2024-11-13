#' Fast Calculation of Matrix Trace and Hilbert Schmidt Norm
#'
#' Provides multithreaded calculation of trace and Hilbert Schmidt Norm of a matrix \eqn{PMP} (where \eqn{P} is a projection matrix and \eqn{M} is real symmetric) without explicitly forming \eqn{PMP}.
#'
#' \eqn{P} here is assumed to have the form \eqn{I-QQ'} for some matrix \eqn{Q} of orthogonal columns.
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
#'        haplotype index at which to start trace calculation --- useful for distributed computation (experimental feature, more documentation to come<!-- TODO -->)
#' @param nthreads
#'        the number of CPU cores to use.
#'        By default uses the `parallel` package to detect the number of physical cores.
#'
#' @return
#' A list containing three elements:
#'
#' \describe{
#'   \item{`trace`}{the trace, \eqn{\mathrm{tr}(PMP)};}
#'   \item{`hsnorm2`}{the *squared* Hilbert Schmidt Norm of \eqn{PMP}, \eqn{\mathrm{tr}((PMP)'PMP)};}
#'   \item{`diag`}{the diagonal of \eqn{PMP}.}
#' }
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
