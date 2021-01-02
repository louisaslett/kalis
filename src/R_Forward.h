#ifndef R_FORWARD_H
#define R_FORWARD_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP ResetForwardTable(SEXP Rfwd);

SEXP Forward(SEXP Rfwd,
             SEXP Rone_step,
             SEXP Rt,
             SEXP RPi,
             SEXP Rmu,
             SEXP Rrho,
             SEXP Ruse_speidel,
             SEXP Rnthreads);

#endif
