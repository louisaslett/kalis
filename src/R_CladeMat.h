#ifndef R_CLADEMAT_H
#define R_CLADEMAT_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP CladeMat(SEXP Rfwd,
              SEXP Rbck,
              SEXP RM,
              SEXP Runit_dist,
              SEXP Rthresh,
              SEXP Rmax1var,
              SEXP Rnthreads);

SEXP UpdateRealInPlace(SEXP RM,
                       SEXP Ridx,
                       SEXP Rvec);

#endif
