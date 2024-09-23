#include "R_CalcTraces.h"

#define _GNU_SOURCE
#include <pthread.h>
#include <math.h>



struct CalcTraces_A_core_args {
  const double *restrict M;
  const double *restrict tX;
  const double *restrict tQ;
  const double *restrict J;
  double *diag;
  size_t r;
  size_t from_off;
  size_t p;
};
struct CalcTraces_A_args {
  struct CalcTraces_A_core_args *core_args;
  double *restrict res;
  double *restrict res2;
  size_t from;
  size_t N;
};

void* CalcTraces_A(void *args) {
  struct CalcTraces_A_args *ct_args;
  ct_args = (struct CalcTraces_A_args *) args;
  const double *restrict M = ct_args->core_args->M;
  const double *restrict tX = ct_args->core_args->tX;
  const double *restrict tQ = ct_args->core_args->tQ;
  const double *restrict J = ct_args->core_args->J;
  double *restrict res = ct_args->res;
  double *restrict res2 = ct_args->res2;
  double *diag = ct_args->core_args->diag;
  size_t r = ct_args->core_args->r;
  size_t from_off = ct_args->core_args->from_off;
  size_t from = ct_args->from;
  size_t N = ct_args->N;
  size_t p = ct_args->core_args->p;

  double temp;
  double temp2;
  for(size_t j = from; j < from+N; j++) {
    for(size_t i = 0; i < r; i++) {
      temp = M[i + j*r];
      for(size_t l = 0; l < p; l++){
        // temp += temp;
        temp += tX[l + i*p] * tQ[l + (from_off+j)*p] - tQ[l + i*p] * J[l + (from_off+j)*p];
      }
      temp2 = temp * temp;
      res2[0] += temp2; // part of the HS norm
      if(i==j){
        res[0] += temp; // add to the trace (for the expected value of the quadratic form)
        diag[j] = temp; // add to the diagonal component of the varaince (for the variance of the quadratic form)
      }
    }
  }

  return(NULL);
}


// // [[Rcpp::export]]
// List CalcTraces(NumericMatrix M,  // r x c
//                 NumericMatrix tX, // p x r
//                 NumericMatrix tQ, // p x r
//                 NumericMatrix J,  // p x r (will only use a subset of these rows)
//                 size_t from_recipient,
//                 size_t nthreads) {
SEXP CalcTraces(SEXP RM,  // r x c
                SEXP RtX, // p x r
                SEXP RtQ, // p x r
                SEXP RJ,  // p x r (will only use a subset of these rows)
                SEXP Rfrom_recipient,
                SEXP Rnthreads) {

  size_t p = (size_t) Rf_nrows(RtX);
  size_t r = (size_t) Rf_nrows(RM);
  size_t c = (size_t) Rf_ncols(RM);

  size_t nthreads = (size_t) Rf_asInteger(Rnthreads);

  SEXP Rtrace = PROTECT(Rf_allocVector(REALSXP, 1));
  REAL(Rtrace)[0] = 0.0;
  SEXP Rhsnorm = PROTECT(Rf_allocVector(REALSXP, 1));
  REAL(Rhsnorm)[0] = 0.0;
  SEXP Rdiag = PROTECT(Rf_allocVector(REALSXP, r));
  for(size_t i=0; i<r; i++) {
    REAL(Rdiag)[i] = 0.0;
  }

  size_t from_off = (size_t) (Rf_asInteger(Rfrom_recipient)-1);

  struct CalcTraces_A_core_args core_args = {
    .M = REAL(RM),
    .tX = REAL(RtX),
    .tQ = REAL(RtQ),
    .J = REAL(RJ),
    .diag = REAL(Rdiag),
    .r = r,
    .from_off = from_off,
    .p = p
  };

  if(nthreads > 1) {
    pthread_t threads[nthreads];
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    double res_perth[nthreads+1]; // = (double*) R_alloc((nthreads+1), sizeof(double));
    double res_perth2[nthreads+1]; // = (double*) R_alloc((nthreads+1), sizeof(double));

    size_t num_perth = c/nthreads;
    size_t rag_end   = c%nthreads;

    struct CalcTraces_A_args args[nthreads+1];
    for(size_t i=0; i<nthreads; i++) {
      args[i].core_args = &core_args;
      args[i].res = res_perth + i;
      args[i].res2 = res_perth2 + i;
      args[i].from = i*num_perth;
      args[i].N = num_perth;
    };

    for(size_t i=0; i<nthreads; ++i) {
      pthread_create(&threads[i], &attr, CalcTraces_A, (void*) &args[i]);
    }
    // Tidy ragged end
    if(rag_end != 0) {
      args[nthreads].core_args = &core_args;
      args[nthreads].res = res_perth + nthreads;
      args[nthreads].res2 = res_perth2 + nthreads;
      args[nthreads].from = nthreads*num_perth;
      args[nthreads].N = rag_end;
      CalcTraces_A((void*) &args[nthreads]);
    }

    for(size_t i=0; i<nthreads; i++) {
      pthread_join(threads[i], NULL);
    }
    pthread_attr_destroy(&attr);

    for(size_t j = 0; j < nthreads+1; j++) {
      REAL(Rtrace)[0] += res_perth[j];
      REAL(Rhsnorm)[0] += res_perth2[j];
    }
  } else {
    struct CalcTraces_A_args args;
    args.core_args = &core_args;
    args.res = REAL(Rtrace);
    args.res2 = REAL(Rhsnorm);
    args.from = 0;
    args.N = c;
    CalcTraces_A((void*) &args);
  }

  const char *names[] = {"trace", "hsnorm2", "diag", ""};
  SEXP RL = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(RL, 0, Rtrace);
  SET_VECTOR_ELT(RL, 1, Rhsnorm);
  SET_VECTOR_ELT(RL, 2, Rdiag);

  UNPROTECT(4);
  return(RL);
}
