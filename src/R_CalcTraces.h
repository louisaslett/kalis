#ifndef R_CALCTRACES_H
#define R_CALCTRACES_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP CalcTraces(SEXP RM,
                SEXP RtX,
                SEXP RtQ,
                SEXP RJ,
                SEXP Rfrom_recipient,
                SEXP Rnthreads);
#endif
