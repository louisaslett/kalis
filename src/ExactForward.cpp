#include <Rcpp.h>
using namespace Rcpp;

#include <stdlib.h>
#include <yepLibrary.h>
#include <yepCore.h>
#include <yepMath.h>

#include "Cache.h"

// [[Rcpp::export]]
NumericMatrix ExactForwardNaiveC_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
  NumericMatrix alpha(N, N);
  double theta;
  char donor_hap, recipient_hap, H;
  int_fast32_t l=0;

  double *f, *fold;
  f    = (double*) malloc(sizeof(double)*N);
  fold = (double*) malloc(sizeof(double)*N);

  // Locus zero setup
  for(int_fast32_t recipient=0; recipient<N; ++recipient) {
    recipient_hap = seq_ind[recipient][0];
    fold[recipient] = 0.0;

    for(int_fast32_t donor=0; donor<N; ++donor) {
      donor_hap = seq_ind[donor][0];
      H = (recipient_hap ^ donor_hap) & 1;
      theta = (H * mu[0]
                 + (1-H) * (1.0 - mu[0]));

      fold[recipient] += alpha(donor, recipient) = theta*Pi(donor, recipient);

      alpha(donor, recipient) = log(alpha(donor, recipient));
    }

    fold[recipient] = -log(fold[recipient]*rho[0]);
  }

  while(l<t) {
    ++l;
    for(int_fast32_t recipient=0; recipient<N; ++recipient) {
      recipient_hap = seq_ind[recipient][l/8];
      f[recipient] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        donor_hap = seq_ind[donor][l/8];
        H = ((recipient_hap ^ donor_hap) >> (l%8)) & 1;
        theta = (H * mu[l]
                   + (1-H) * (1.0 - mu[l]));

        f[recipient] += alpha(donor, recipient) = (theta * Pi(donor, recipient)
                                                     + theta * (1.0-rho[l-1]) * exp(alpha(donor, recipient) + fold[recipient]));

        alpha(donor, recipient) = log(alpha(donor, recipient)) - fold[recipient];
      }

      fold[recipient] = -(log(f[recipient] * rho[l]) - fold[recipient]);
    }
  }

  free(f);
  free(fold);

  return(alpha);
}

#include "ExactForwardISPC.h"
// [[Rcpp::export]]
NumericMatrix dExactForward_ISPC_st_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
  NumericMatrix alpha(N, N);
  // float *alphaf = (float*) malloc(sizeof(float)*N*N),
  //   *Pif = (float*) malloc(sizeof(float)*N*N),
  //   *muf = (float*) malloc(sizeof(float)*N),
  //   *rhof = (float*) malloc(sizeof(float)*N);
  //
  // for(int i=0; i<N; i++) {
  //   muf[i] = (float) mu(i);
  //   rhof[i] = (float) rho(i);
  //   for(int j=0; j<N; j++) {
  //     Pif[i*N+j] = (float) Pi(j,i);
  //   }
  // }

  // ispc::dExactForward_ISPC_st(t, L, N, Pif, muf, rhof, alphaf);
  ispc::dExactForward_ISPC_st(t, L, N, &(Pi(0,0)), &(mu(0)), &(rho(0)), &(alpha(0,0)));

  // for(int j=0; j<N; j++) {
  //   for(int i=0; i<N; i++) {
  //     alpha(i,j) = (double) alphaf[i+j*N];
  //   }
  // }

  return(alpha);
}

// [[Rcpp::export]]
NumericMatrix ExactForwardYepppExpC_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
  NumericMatrix alpha(N, N);
  double theta;
  char donor_hap, recipient_hap, H;
  int_fast32_t l=0;

  yepLibrary_Init();

  double *f, *fold;
  f    = (double*) malloc(sizeof(double)*N);
  fold = (double*) malloc(sizeof(double)*N);

  // Locus zero setup
  for(int_fast32_t recipient=0; recipient<N; ++recipient) {
    recipient_hap = seq_ind[recipient][0];
    fold[recipient] = 0.0;

    for(int_fast32_t donor=0; donor<N; ++donor) {
      donor_hap = seq_ind[donor][0];
      H = (recipient_hap ^ donor_hap) & 1;
      theta = (H * mu[0]
                 + (1-H) * (1.0 - mu[0]));

      fold[recipient] += alpha(donor, recipient) = theta*Pi(donor, recipient);

      alpha(donor, recipient) = log(alpha(donor, recipient));
    }

    fold[recipient] = -log(fold[recipient]*rho[0]);
  }

  while(l<t) {
    ++l;
    for(int_fast32_t recipient=0; recipient<N; ++recipient) {
      recipient_hap = seq_ind[recipient][l/8];
      f[recipient] = 0.0;

      double *alphaRow = &(alpha[N*recipient]);
      double YepppTmp[N];
      yepCore_Add_V64fS64f_V64f(&(alpha(0, recipient)), fold[recipient], YepppTmp, N);
      yepMath_Exp_V64f_V64f(YepppTmp, alphaRow, N);

      double *PiRow = &(Pi(0, recipient));
      double muTmp = mu[l];
      double *rhoTmp = &(rho[l-1]);
      for(int_fast32_t donor=0; donor<N; ++donor) {
        donor_hap = seq_ind[donor][l/8];
        H = ((recipient_hap ^ donor_hap) >> (l%8)) & 1;
        theta = (H * muTmp
                   + (1-H) * (1.0 - muTmp));

        f[recipient] += alphaRow[donor] = (theta * PiRow[donor]
                                                     + theta * (1.0-rhoTmp[0]) * alphaRow[donor]);
      }
      yepMath_Log_V64f_V64f(alphaRow, YepppTmp, N);
      yepCore_Subtract_V64fS64f_V64f(YepppTmp, fold[recipient], alphaRow, N);

      fold[recipient] = -(log(f[recipient] * rhoTmp[1]) - fold[recipient]);
    }
  }

  free(f);
  free(fold);

  yepLibrary_Release();

  return(alpha);
}
