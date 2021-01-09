#ifndef R_MATANDMUL_H
#define R_MATANDMUL_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP MatAndMul(SEXP RM,
               SEXP Rfwd,
               SEXP Rbck,
               SEXP Rx,
               SEXP Rstd,
               SEXP Rcalc_probs,
               SEXP Runif_underflow,
               SEXP Rnthreads);

#endif
