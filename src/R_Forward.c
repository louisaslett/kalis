#include "R_Forward.h"

#include "R_Kalis.h"
#include "Stencil.h"
#include "Cache.h"
#include "ExactForward.h"

SEXP ResetForwardTable(SEXP Rfwd) {
  int *l;

  KALIS_GET_TABLE_POS(l, Rfwd);

  *l = 2147483647;

  KALIS_RETURN;
}

SEXP Forward(SEXP Rfwd,
             SEXP Rone_step,
             SEXP Rt,
             SEXP RPi,
             SEXP Rmu,
             SEXP Rrho,
             SEXP Ruse_speidel,
             SEXP Rnthreads) {
  const int L = hap_size;
  const int N = num_inds;
  double *alpha, *alpha_f, *Pi, *mu, *rho;
  int *one_step, *l, *from_rec, *to_rec, *t, *use_speidel, *nthreads, nthreads_n;
  int which_mu, which_Pi, which_speidel, which_step;

  // Extract forward table variables
  KALIS_GET_TABLE(alpha, Rfwd);
  KALIS_GET_TABLE_AUX(alpha_f, Rfwd);
  KALIS_GET_TABLE_POS(l, Rfwd);
  KALIS_GET_TABLE_FROM(from_rec, Rfwd);
  KALIS_GET_TABLE_TO(to_rec, Rfwd);

  // Extract other function arguments
  one_step = LOGICAL(Rone_step);
  if(XLENGTH(Rone_step) != 1) {
    REprintf("Error: one_step should be scalar.\n");
    KALIS_RETURN;
  }
  if(*one_step) {
    which_step = ONESTEP;
  } else {
    which_step = NSTEP;
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
  if(*l<=L && *l>*t) {
    REprintf("Error: The forward table provided is for locus position %d which is already past requested locus %d.\n", *l, *t);
    KALIS_RETURN;
  }
  if(*l==*t) {
    KALIS_RETURN;
  }

  // Dispatch to correct optimised core function
  switch(which_mu + which_Pi + which_speidel + which_step) {
  case MU_SCALAR + PI_SCALAR + NOSPEIDEL + NSTEP:
    ExactForward_scalarMu_scalarPi(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_SCALAR + NOSPEIDEL + NSTEP:
    ExactForward_scalarPi(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, mu, rho, nthreads, nthreads_n); break;
  case MU_SCALAR + PI_MATRIX + NOSPEIDEL + NSTEP:
    ExactForward_scalarMu(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_MATRIX + NOSPEIDEL + NSTEP:
    ExactForward(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, Pi, mu, rho, nthreads, nthreads_n); break;
  case MU_SCALAR + PI_SCALAR + SPEIDEL + NSTEP:
    ExactForward_speidel_scalarMu_scalarPi(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_SCALAR + SPEIDEL + NSTEP:
    ExactForward_speidel_scalarPi(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, mu, rho, nthreads, nthreads_n); break;
  case MU_SCALAR + PI_MATRIX + SPEIDEL + NSTEP:
    ExactForward_speidel_scalarMu(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_MATRIX + SPEIDEL + NSTEP:
    ExactForward_speidel(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, Pi, mu, rho, nthreads, nthreads_n); break;
  case MU_SCALAR + PI_SCALAR + NOSPEIDEL + ONESTEP:
    ExactForward_1step_scalarMu_scalarPi(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_SCALAR + NOSPEIDEL + ONESTEP:
    REprintf("Error: 1-step forward currently only available for scalar mu and scalar Pi.\n"); KALIS_RETURN;
  case MU_SCALAR + PI_MATRIX + NOSPEIDEL + ONESTEP:
    REprintf("Error: 1-step forward currently only available for scalar mu and scalar Pi.\n"); KALIS_RETURN;
  case MU_VECTOR + PI_MATRIX + NOSPEIDEL + ONESTEP:
    REprintf("Error: 1-step forward currently only available for scalar mu and scalar Pi.\n"); KALIS_RETURN;
  case MU_SCALAR + PI_SCALAR + SPEIDEL + ONESTEP:
    ExactForward_1step_speidel_scalarMu_scalarPi(
      alpha, alpha_f, *from_rec-1, *l-1, *t-1,
      *from_rec-1, *to_rec, L, N, *Pi, *mu, rho, nthreads, nthreads_n); break;
  case MU_VECTOR + PI_SCALAR + SPEIDEL + ONESTEP:
    REprintf("Error: 1-step forward currently only available for scalar mu and scalar Pi.\n"); KALIS_RETURN;
  case MU_SCALAR + PI_MATRIX + SPEIDEL + ONESTEP:
    REprintf("Error: 1-step forward currently only available for scalar mu and scalar Pi.\n"); KALIS_RETURN;
  case MU_VECTOR + PI_MATRIX + SPEIDEL + ONESTEP:
    REprintf("Error: 1-step forward currently only available for scalar mu and scalar Pi.\n"); KALIS_RETURN;
  default:
    REprintf("Error: unknown combination of runtime options.\n"); KALIS_RETURN;
  }

  *l = *t;

  KALIS_RETURN;
}
