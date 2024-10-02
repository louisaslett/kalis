#ifndef R_OPTCKPT_H
#define R_OPTCKPT_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP OptCkpt(SEXP Rcost_table, SEXP Rindex_table, SEXP Rpropagation_cost);

#endif
