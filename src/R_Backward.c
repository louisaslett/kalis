#include "R_Backward.h"

#include "R_Kalis.h"
#include "Stencil.h"
#include "Cache.h"
#include "ExactBackward.h"

SEXP ResetBackwardTable(SEXP Rbck) {
  int *l, *beta_theta;

  KALIS_GET_TABLE_POS(l, Rbck);
  KALIS_GET_TABLE_BETATHETA(beta_theta, Rbck);

  *l = 2147483647;
  *beta_theta = 0;

  KALIS_RETURN;
}

SEXP Backward(SEXP Rbck,
              SEXP Rend_beta_theta,
              SEXP Rt,
              SEXP RPi,
              SEXP Rmu,
              SEXP Rrho,
              SEXP Ruse_speidel,
              SEXP Rnthreads) {
  const int L = hap_size;
  const int N = num_inds;
  double *beta, *beta_g, *Pi, *mu, *rho;
  int *cur_beta_theta, *end_beta_theta, *l, *from_rec, *to_rec, *t, *use_speidel, *nthreads, nthreads_n;
  int which_mu, which_Pi, which_speidel;

  // Extract backward table variables
  KALIS_GET_TABLE(beta, Rbck);
  KALIS_GET_TABLE_AUX(beta_g, Rbck);
  KALIS_GET_TABLE_BETATHETA(cur_beta_theta, Rbck);
  KALIS_GET_TABLE_POS(l, Rbck);
  KALIS_GET_TABLE_FROM(from_rec, Rbck);
  KALIS_GET_TABLE_TO(to_rec, Rbck);

  // Extract other function arguments
  end_beta_theta = LOGICAL(Rend_beta_theta);
  if(XLENGTH(Rend_beta_theta) != 1) {
    REprintf("Error: end_beta_theta should be scalar.\n");
    KALIS_RETURN;
  }
  t = INTEGER(Rt);
  if(XLENGTH(Rt) != 1) {
    REprintf("Error: t should be scalar.\n");
    KALIS_RETURN;
  }
  Pi = REAL(RPi);
  if(Rf_isMatrix(RPi) && Rf_nrows(RPi) == N && Rf_ncols(RPi) == N) {
    which_Pi = PI_MATRIX;
  } else if(!Rf_isMatrix(RPi) && XLENGTH(RPi) == 1) {
    which_Pi = PI_SCALAR;
  } else {
    REprintf("Error: Pi should be N x N matrix or scalar.\n");
    KALIS_RETURN;
  }
  mu = REAL(Rmu);
  if(XLENGTH(Rmu) == L) {
    which_mu = MU_VECTOR;
  } else if(XLENGTH(Rmu) == 1) {
    which_mu = MU_SCALAR;
  } else {
    REprintf("Error: mu should be length L or scalar.\n");
    KALIS_RETURN;
  }
  rho = REAL(Rrho);
  if(XLENGTH(Rrho) != L) {
    REprintf("Error: rho should be length L.\n");
    KALIS_RETURN;
  }
  use_speidel = LOGICAL(Ruse_speidel);
  if(XLENGTH(Ruse_speidel) != 1) {
    REprintf("Error: use_speidel should be scalar.\n");
    KALIS_RETURN;
  }
  if(*use_speidel) {
    which_speidel = SPEIDEL;
  } else {
    which_speidel = NOSPEIDEL;
  }
  nthreads = INTEGER(Rnthreads);
  nthreads_n = LENGTH(Rnthreads);
  if(nthreads_n > 1) {
#if !defined(KALIS_AFFINITY)
    Rprintf("Thread affinity not supported on this platform, running on %d cores without setting affinity.\n", nthreads_n);
#endif
  }

  // Check sanity of requested locus
  if(*l<*t) {
    REprintf("Error: The backward table provided is for locus position %d which is already past requested locus %d.\n", *l, *t);
    KALIS_RETURN;
  }
  if(*l==*t && !(!*cur_beta_theta && *end_beta_theta)) {
    KALIS_RETURN;
  }

  // Dispatch to correct optimised core function
  switch(which_mu + which_Pi + which_speidel) {
  case MU_SCALAR + PI_SCALAR + NOSPEIDEL:
    ExactBackward_scalarMu_scalarPi(
      beta, beta_g, cur_beta_theta, end_beta_theta, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_SCALAR + NOSPEIDEL:
    ExactBackward_scalarPi(
      beta, beta_g, cur_beta_theta, end_beta_theta, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, mu, rho, nthreads, nthreads_n); break;
  case MU_SCALAR + PI_MATRIX + NOSPEIDEL:
    ExactBackward_scalarMu(
      beta, beta_g, cur_beta_theta, end_beta_theta, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_MATRIX + NOSPEIDEL:
    ExactBackward(
      beta, beta_g, cur_beta_theta, end_beta_theta, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, Pi, mu, rho, nthreads, nthreads_n); break;
  case MU_SCALAR + PI_SCALAR + SPEIDEL:
    ExactBackward_speidel_scalarMu_scalarPi(
      beta, beta_g, cur_beta_theta, end_beta_theta, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_SCALAR + SPEIDEL:
    ExactBackward_speidel_scalarPi(
      beta, beta_g, cur_beta_theta, end_beta_theta, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, mu, rho, nthreads, nthreads_n); break;
  case MU_SCALAR + PI_MATRIX + SPEIDEL:
    ExactBackward_speidel_scalarMu(
      beta, beta_g, cur_beta_theta, end_beta_theta, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_MATRIX + SPEIDEL:
    ExactBackward_speidel(
      beta, beta_g, cur_beta_theta, end_beta_theta, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, Pi, mu, rho, nthreads, nthreads_n); break;
  default:
    REprintf("Error: unknown combination of runtime options.\n"); KALIS_RETURN;
  }

  *l = *t;

  KALIS_RETURN;
}
