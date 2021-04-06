#include "R_MatAndMulBtwVar.h"

#define _GNU_SOURCE
#include <pthread.h>
#include <math.h>

#include "R_Kalis.h"



void MatAndMulBtwVar_C_probs(const double* restrict alpha,
                             const double* restrict f,
                             const double* restrict beta,
                             const double* restrict g,
                             const double* restrict x,
                             double* restrict c_1,
                             double* restrict res,
                             int useunif,
                             size_t j,
                             size_t r,
                             size_t from_off,
                             double fwd_rho,
                             double bck_rho) {
  double z0 = 0.0;
  double rd = (double) (r - 1);
  double ccc = 0.0;

  if(useunif) {
    ccc = 1/rd;
  }

  if((*g != 0) & (*f != 0)) {
    double beta_coeff_c1 = (1 - bck_rho) / (*g);
    double alpha_coeff_c1 = (1 - fwd_rho) / (*f);

    double alpha_offset = fwd_rho / rd;

    // We will assume in this accumulation that alpha * beta is always 0 along the matrix diagonal
    for(size_t i = 0; i < r; i++) {
      z0 += c_1[i] = (alpha_coeff_c1 * *(alpha++) + alpha_offset) * (beta_coeff_c1 * *(beta++) + bck_rho) ;
    }
  }

  if(z0 <= 0) {
    // if the sum is zero (all of the alpha*beta entries = 0), then all of the distances are beyond the resolution
    // of numerical precision, so we set all of the distances to the -log of the smallest respresentable double
    for(size_t i = 0; i < r; i++) {
      if(i==j + from_off) {
        c_1[i] = 0.0;
      } else {
        c_1[i] = ccc;
      }
    }

  } else { // z0 > 0

    for(size_t i = 0; i < r; i++) {
      if(i==j + from_off) {
        c_1[i] = 0.0;
      } else {

        c_1[i] = c_1[i] / z0;

        res[j+from_off] += c_1[i]*x[i];
        res[i]          += c_1[i]*x[j+from_off];
      }
    }
  }
}

void MatAndMulBtwVar_C_raw_dist(const double* restrict alpha,
                                const double* restrict f,
                                const double* restrict beta,
                                const double* restrict g,
                                const double* restrict x,
                                double* restrict c_1,
                                double* restrict res,
                                size_t j,
                                size_t r,
                                size_t from_off,
                                double fwd_rho,
                                double bck_rho,
                                int minus_min) {
  double z0 = 0.0;
  double rd = (double) (r - 1);

  if((*g != 0) & (*f != 0)) {
    double beta_coeff_c1 = (1 - bck_rho) / (*g);
    double alpha_coeff_c1 = (1 - fwd_rho) / (*f);

    double alpha_offset = fwd_rho / rd;

    // We will assume in this accumulation that alpha * beta is always 0 along the matrix diagonal
    if(minus_min) {
      for(size_t i = 0; i < r; i++) {
        c_1[i] = (alpha_coeff_c1 * *(alpha++) + alpha_offset) * (beta_coeff_c1 * *(beta++) + bck_rho);
        z0 = (z0>c_1[i])?z0:c_1[i];
      }
    } else {
      for(size_t i = 0; i < r; i++) {
        z0 += c_1[i] = (alpha_coeff_c1 * *(alpha++) + alpha_offset) * (beta_coeff_c1 * *(beta++) + bck_rho) ;
      }
    }
  }

  if(z0 <= 0) {
    // if the sum is zero (all of the alpha*beta entries = 0), then all of the distances are beyond the resolution
    // of numerical precision, so we set all of the distances to the -log of the smallest respresentable double
    for(size_t i = 0; i < r; i++) {
      if(i==j + from_off) {
        c_1[i] = 0.0;
      } else {
        c_1[i] = 744.4400719213812180897;
      }
    }

  } else { // z0 > 0

    for(size_t i = 0; i < r; i++) {
      if(i==j + from_off) {
        c_1[i] = 0.0;
      } else {

        c_1[i] = c_1[i] / z0;

        if(c_1[i] == 0) {
          c_1[i] = 744.4400719213812180897;  // -log of smallest representable double
        } else {
          c_1[i] = -log(c_1[i]); // this is the final value for M[i,j]
        }
        res[j+from_off] += c_1[i]*x[i];
        res[i]          += c_1[i]*x[j+from_off];
      }
    }

  }
}



void MatAndMulBtwVar_C_standardize_dist(const double* restrict alpha,
                                        const double* restrict f,
                                        const double* restrict beta,
                                        const double* restrict g,
                                        const double* restrict x,
                                        double* restrict c_1,
                                        double* restrict res,
                                        size_t j,
                                        size_t r,
                                        size_t from_off,
                                        double fwd_rho,
                                        double bck_rho) {
  double z0 = 0.0;
  double z1 = 0.0;
  double z2 = 0.0;
  double rd = (double) (r - 1);

  if((*g != 0) & (*f != 0)) {

    double beta_coeff_c1 = (1 - bck_rho) / (*g);
    double alpha_coeff_c1 = (1 - fwd_rho) / (*f);
    double alpha_offset = fwd_rho / rd;


    // We will assume in this accumulation that alpha * beta is always 0 along the matrix diagonal
    for(size_t i = 0; i < r; i++) {
      z0 += c_1[i] = (alpha_coeff_c1 * *(alpha++) + alpha_offset) * (beta_coeff_c1 * *(beta++) + bck_rho) ;
    }
  }

  // if all of the donors have "0" probability of being copied, so all c_1[i] = 0, we exit here...otherwise...
  if(z0<=0) {

    for(size_t i = 0; i < r; i++) {
      c_1[i] = 0.0;
    }

  } else { // z0 > 0

    for(size_t i=0; i < r; i++) {
      if(i==j + from_off) {
        continue;
      } else {
        c_1[i] = c_1[i] / z0;

        if(c_1[i]==0) {
          c_1[i] = 744.4400719213812180897;  // -log of smallest representable double
        } else {
          c_1[i] = -log(c_1[i]);
        }

        z1 += c_1[i];
        z2 += c_1[i]*c_1[i];
      }
    }

    z1 = z1 / rd; // Calculate Mean

    z2 = (z2 - z1*z1 * rd) / (rd - 1.0); // Unbiased Estimate of the Variance

    if(z2 <= 0) { // if there is no variance in the distances, we just set them all to 0
      for(size_t i = 0; i < r; i++) {
        c_1[i] = 0.0;
      }
    } else {

      z2 = 1/sqrt(z2); // Calculate standard deviation

      z1 = -z1 * z2;

      for(size_t i = 0; i < r; i++) {
        if(i==j + from_off) {
          c_1[i] = 0.0;
        } else {
          c_1[i] = c_1[i] * z2 + z1; // this is the final value for M[i,j]
          res[j+from_off] += c_1[i]*x[i];
          res[i]          += c_1[i]*x[j+from_off];
        }
      }
    }

  }
}



struct MatAndMulBtwVar_B_core_args {
  double* M;
  const double* restrict alpha;
  const double* restrict f;
  const double* restrict beta;
  const double* restrict g;
  const double* restrict x;
  int standardize;
  int minus_min;
  int probs;
  int useunif;
  size_t r;
  size_t from_off;
  double fwd_rho;
  double bck_rho;
};
struct MatAndMulBtwVar_B_args {
  struct MatAndMulBtwVar_B_core_args *core_args;
  double* restrict res;
  size_t from;
  size_t N;
};

void* MatAndMulBtwVar_B(void *args) {
  struct MatAndMulBtwVar_B_args *mat_args;
  mat_args = (struct MatAndMulBtwVar_B_args *) args;
  double* M = mat_args->core_args->M;
  const double* restrict alpha = mat_args->core_args->alpha;
  const double* restrict f = mat_args->core_args->f;
  const double* restrict beta = mat_args->core_args->beta;
  const double* restrict g = mat_args->core_args->g;
  const double* restrict x = mat_args->core_args->x;
  int standardize = mat_args->core_args->standardize;
  int minus_min = mat_args->core_args->minus_min;
  int probs = mat_args->core_args->probs;
  int useunif = mat_args->core_args->useunif;
  size_t r = mat_args->core_args->r;
  size_t from_off = mat_args->core_args->from_off;
  double fwd_rho = mat_args->core_args->fwd_rho;
  double bck_rho = mat_args->core_args->bck_rho;
  double* restrict res = mat_args->res;
  size_t from = mat_args->from;
  size_t N = mat_args->N;

  if(probs) {

    for(size_t j = from; j < from+N; j++) {
      MatAndMulBtwVar_C_probs(alpha+j*r,
                              f+j,
                              beta+j*r,
                              g+j,
                              x,
                              M+j*r,
                              res,
                              useunif,
                              j,
                              r,
                              from_off,
                              fwd_rho,
                              bck_rho);
    }

  } else {

    if(standardize) {
      for(size_t j = from; j < from+N; j++) {
        MatAndMulBtwVar_C_standardize_dist(alpha+j*r,
                                           f+j,
                                           beta+j*r,
                                           g+j,
                                           x,
                                           M+j*r,
                                           res,
                                           j,
                                           r,
                                           from_off,
                                           fwd_rho,
                                           bck_rho);
      }
    } else {
      for(size_t j = from; j < from+N; j++) {
        MatAndMulBtwVar_C_raw_dist(alpha+j*r,
                                   f+j,
                                   beta+j*r,
                                   g+j,
                                   x,
                                   M+j*r,
                                   res,
                                   j,
                                   r,
                                   from_off,
                                   fwd_rho,
                                   bck_rho,
                                   minus_min);
      }
    }

  }

  return(NULL);
}



void MatAndMulBtwVar_A(double* restrict res,
                       double* restrict M,
                       const double* restrict alpha,
                       const double* restrict f,
                       const double* restrict beta,
                       const double* restrict g,
                       const double* restrict x,
                       int standardize,
                       int minus_min,
                       int probs,
                       int useunif,
                       size_t from_off,
                       size_t nthreads,
                       size_t r,
                       size_t c,
                       double fwd_rho,
                       double bck_rho) {

  from_off = (from_off-1); // got rid of /2 here

  struct MatAndMulBtwVar_B_core_args core_args = {
    .M = M,
    .alpha = alpha,
    .f = f,
    .beta = beta,
    .g = g,
    .x = x,
    .standardize = standardize,
    .minus_min = minus_min,
    .probs = probs,
    .useunif = useunif,
    .r = r,
    .from_off = from_off,
    .fwd_rho = fwd_rho,
    .bck_rho = bck_rho
  };

  if(nthreads > 1) {

    pthread_t threads[nthreads];
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    double *res_perth = (double*) R_alloc(r*(nthreads+1), sizeof(double));
    if (res_perth == NULL) {
      printf("Failed allocating res_perth!\n");
      exit(1);
    }

    size_t num_perth = c/nthreads;
    size_t rag_end   = c%nthreads;

    struct MatAndMulBtwVar_B_args args[nthreads+1];
    for(size_t i=0; i<nthreads; i++) {
      args[i].core_args = &core_args;
      args[i].res = res_perth + i*r;
      args[i].from = i*num_perth;
      args[i].N = num_perth;
    }

    for(size_t i=0; i<nthreads; ++i) {
      pthread_create(&threads[i], &attr, MatAndMulBtwVar_B, (void*) &args[i]);
    }
    // Tidy ragged end
    if(rag_end != 0) {
      args[nthreads].core_args = &core_args;
      args[nthreads].res = res_perth+nthreads*r;
      args[nthreads].from = nthreads*num_perth;
      args[nthreads].N = rag_end;
      MatAndMulBtwVar_B((void*) &args[nthreads]);
    }

    for(size_t i=0; i<nthreads; i++) {
      pthread_join(threads[i], NULL);
    }
    pthread_attr_destroy(&attr);

    for(size_t i = 0; i < r; i++) {
      for(size_t j = 0; j < nthreads+1; j++) {
        res[i] += res_perth[i+j*r];
      }
      res[i] *= 0.5;
    }

  } else {

    struct MatAndMulBtwVar_B_args args;
    args.core_args = &core_args;
    args.res = res;
    args.from = 0;
    args.N = c;

    MatAndMulBtwVar_B((void*) &args);

    for(size_t i = 0; i < r; i++) {
      res[i] *= 0.5;
    }

  }
}



SEXP MatAndMulBtwVar(SEXP RM,
                     SEXP Rfwd,
                     SEXP Rbck,
                     SEXP Rx,
                     SEXP Rstdz,
                     SEXP Rminus_min,
                     SEXP Rcalc_probs,
                     SEXP Runif_underflow,
                     SEXP Rfwd_rho,
                     SEXP Rbck_rho,
                     SEXP Rnthreads) {
  double *restrict alpha, *restrict f, *restrict beta, *restrict g, *restrict M, *restrict x, *restrict res;
  int *restrict from_rec;

  KALIS_GET_TABLE(alpha, Rfwd);
  KALIS_GET_TABLE_AUX(f, Rfwd);
  KALIS_GET_TABLE(beta, Rbck);
  KALIS_GET_TABLE_AUX(g, Rbck);
  KALIS_GET_TABLE_FROM(from_rec, Rfwd);

  size_t r = (size_t) Rf_nrows(VECTOR_ELT(Rfwd, KALIS_IDX_TABLE)); // we need to make p=r and c2 = c
  size_t c = (size_t) Rf_ncols(VECTOR_ELT(Rfwd, KALIS_IDX_TABLE));

  if(Rf_ncols(VECTOR_ELT(Rfwd, KALIS_IDX_TABLE)) != Rf_ncols(VECTOR_ELT(Rbck, KALIS_IDX_TABLE)) ||
     Rf_nrows(VECTOR_ELT(Rfwd, KALIS_IDX_TABLE)) != Rf_nrows(VECTOR_ELT(Rbck, KALIS_IDX_TABLE))) {
    REprintf("Error: alpha, beta aren't right!\n");
    KALIS_RETURN;
  }

  if((size_t) XLENGTH(Rx) != r) {
    REprintf("Error: x isn't right!\n");
    KALIS_RETURN;
  }
  x = REAL(Rx);

  if(!Rf_isMatrix(RM) || (size_t) Rf_nrows(RM) != r || (size_t) Rf_ncols(RM) != c) {
    REprintf("Error: M isn't right!\n");
    KALIS_RETURN;
  }
  M = REAL(RM);

  SEXP Mattrib = PROTECT(Rf_allocVector(STRSXP, 2));
  SEXP Mclass1 = PROTECT(Rf_mkChar("kalisDistanceMatrix"));
  SEXP Mclass2 = PROTECT(Rf_mkChar("matrix"));
  SET_STRING_ELT(Mattrib, 0, Mclass1);
  SET_STRING_ELT(Mattrib, 1, Mclass2);
  Rf_setAttrib(RM, R_ClassSymbol, Mattrib);

  int stdz = Rf_asLogical(Rstdz);
  int minus_min = Rf_asLogical(Rminus_min);
  int probs = Rf_asLogical(Rcalc_probs);
  int useunif = Rf_asLogical(Runif_underflow);

  if(stdz && r < 3) {
    REprintf("Error: each matrix must have at least three rows to perform standardization\n");
    KALIS_RETURN;
  }

  SEXP Rres = PROTECT(Rf_allocVector(REALSXP, r));
  res = REAL(Rres);
  for(size_t i=0; i<r; i++) {
    res[i] = 0.0;
  }

  MatAndMulBtwVar_A(res,
                    M,
                    alpha,
                    f,
                    beta,
                    g,
                    x,
                    stdz,
                    minus_min,
                    probs,
                    useunif,
                    (size_t) *from_rec,
                    (size_t) Rf_asInteger(Rnthreads),
                    r,
                    c,
                    Rf_asReal(Rfwd_rho),
                    Rf_asReal(Rbck_rho));

  UNPROTECT(4);
  return(Rres);
}
