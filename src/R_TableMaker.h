#ifndef R_TABLEMAKER_H
#define R_TABLEMAKER_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP MakeForwardTable(SEXP Rfrom_r, SEXP Rto_r, SEXP Rpars_sha);

SEXP MakeBackwardTable(SEXP Rfrom_r, SEXP Rto_r, SEXP Rpars_sha);

#endif
