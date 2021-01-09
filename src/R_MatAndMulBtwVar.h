#ifndef R_MATANDMULBTWVAR_H
#define R_MATANDMULBTWVAR_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP MatAndMulBtwVar(SEXP RM,
                     SEXP Rfwd,
                     SEXP Rbck,
                     SEXP Rx,
                     SEXP Rstdz,
                     SEXP Rminus_min,
                     SEXP Rcalc_probs,
                     SEXP Runif_underflow,
                     SEXP Rfwd_rho,
                     SEXP Rbck_rho,
                     SEXP Rnthreads);
#endif
