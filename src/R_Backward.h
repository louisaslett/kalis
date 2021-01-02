#ifndef R_BACKWARD_H
#define R_BACKWARD_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP ResetBackwardTable(SEXP Rbck);

SEXP Backward(SEXP Rbck,
              SEXP Rend_beta_theta,
              SEXP Rt,
              SEXP RPi,
              SEXP Rmu,
              SEXP Rrho,
              SEXP Ruse_speidel,
              SEXP Rnthreads);

#endif
