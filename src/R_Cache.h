#ifndef R_CACHE_H
#define R_CACHE_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP ClearHaplotypeCache2();

SEXP CacheHaplotypes_matrix_2(SEXP Rx, SEXP RN, SEXP RL, SEXP Rtranspose);

SEXP CacheHaplotypes_hdf5_2(SEXP Rnexthaps, SEXP Rnexthapsenv, SEXP RN, SEXP RL);

SEXP CacheHaplotypes_hapgz_ncols(SEXP Rfile);

SEXP CacheHaplotypes_hapgz_nlines(SEXP Rfile);

SEXP CacheHaplotypes_hapgz_2(SEXP Rfile, SEXP Rloci_idx, SEXP Rhap_idx, SEXP RL, SEXP RN);

SEXP QueryCache2_ind(SEXP Ridx);

SEXP QueryCache2_loc(SEXP Ridx);

#endif
