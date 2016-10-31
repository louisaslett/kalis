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
  int_fast32_t l=0;

  double *f, *fold;
  f    = (double*) malloc(sizeof(double)*N);
  fold = (double*) malloc(sizeof(double)*N);

  // Locus zero setup
  for(int_fast32_t recipient=0; recipient<N; ++recipient) {
    int32_t recipient_hap = (seq_locus[0][recipient/32] >> recipient%32) & 1;
    fold[recipient] = 0.0;

    for(int_fast32_t donor=0; donor<N; ++donor) {
      int32_t donor_hap = (seq_locus[0][donor/32] >> donor%32) & 1;
      int32_t H = (recipient_hap ^ donor_hap) & 1;
      double theta = (H * mu[0]
                        + (1-H) * (1.0 - mu[0]));

      fold[recipient] += alpha(donor, recipient) = theta*Pi(donor, recipient);

      alpha(donor, recipient) = log(alpha(donor, recipient));
    }

    fold[recipient] = -log(fold[recipient]*rho[0]);
  }

  while(l<t) {
    ++l;
    for(int_fast32_t recipient=0; recipient<N; ++recipient) {
      int32_t recipient_hap = (seq_locus[l][recipient/32] >> recipient%32) & 1;
      f[recipient] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (seq_locus[l][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu[l]
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
  int_fast32_t l=0;

  yepLibrary_Init();

  double *f, *fold;
  f    = (double*) malloc(sizeof(double)*N);
  fold = (double*) malloc(sizeof(double)*N);

  // Locus zero setup
  for(int_fast32_t recipient=0; recipient<N; ++recipient) {
    int32_t recipient_hap = (seq_locus[0][recipient/32] >> recipient%32) & 1;
    fold[recipient] = 0.0;

    for(int_fast32_t donor=0; donor<N; ++donor) {
      int32_t donor_hap = (seq_locus[0][donor/32] >> donor%32) & 1;
      int32_t H = (recipient_hap ^ donor_hap) & 1;
      double theta = (H * mu[0]
                        + (1-H) * (1.0 - mu[0]));

      fold[recipient] += alpha(donor, recipient) = theta*Pi(donor, recipient);

      alpha(donor, recipient) = log(alpha(donor, recipient));
    }

    fold[recipient] = -log(fold[recipient]*rho[0]);
  }

  while(l<t) {
    ++l;
    for(int_fast32_t recipient=0; recipient<N; ++recipient) {
      int32_t recipient_hap = 0;
      recipient_hap -= (seq_locus[l][recipient/32] >> recipient%32) & 1;

      f[recipient] = 0.0;

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *alphaRow = &(alpha[N*recipient]);
      double YepppTmp[N];
      yepCore_Add_V64fS64f_V64f(alphaRow, fold[recipient], YepppTmp, N);
      yepMath_Exp_V64f_V64f(YepppTmp, alphaRow, N);

      // theta[donor] * (1-rho) * alphaRow[donor]
      yepCore_Multiply_IV64fS64f_IV64f(alphaRow, (1.0-rho[l-1]), N);
      double muTmp = mu[l];
      for(int_fast32_t donoroff=0; donoroff<N/32; ++donoroff) {
        int32_t H = (recipient_hap ^ seq_locus[l][donoroff]); // _mm256_xor_si256
        for(char donor=0; donor<32; ++donor) {
          // donor_hap = (seq_locus[l][donoroff] >> donor) & 1;
          // H = (recipient_hap ^ donor_hap) & 1;
          // YepppTmp[donoroff*8+donor] = (H * muTmp
          //                      + (1-H) * (1.0 - muTmp));
          YepppTmp[donoroff*32+donor] = (H&1) * (2 * muTmp - 1.0) - muTmp + 1.0;
          H >>= 1;
        }
      }
      for(char donor=0; donor<N%32; ++donor) {
        int32_t donor_hap = (seq_locus[l][N/32] >> donor) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        YepppTmp[(N/32)*32+donor] = (H * muTmp
                                       + (1-H) * (1.0 - muTmp));
      }
      yepCore_Multiply_IV64fV64f_IV64f(alphaRow, YepppTmp, N);

      // TODO: benchmark against a fused multiply and add for these two together
      // theta[donor] * Pi(donor, recipient)
      yepCore_Multiply_IV64fV64f_IV64f(YepppTmp, &(Pi(0, recipient)), N);
      // Tot up alphaRow
      yepCore_Add_IV64fV64f_IV64f(alphaRow, YepppTmp, N);

      yepCore_Sum_V64f_S64f(alphaRow, f+recipient, N);
      // f[recipient] += alphaRow[donor] = (theta * PiRow[donor]
      //                                      + theta * (1.0-rhoTmp[0]) * alphaRow[donor]);

      yepMath_Log_V64f_V64f(alphaRow, YepppTmp, N);
      yepCore_Subtract_V64fS64f_V64f(YepppTmp, fold[recipient], alphaRow, N);

      fold[recipient] = -(log(f[recipient] * rho[l]) - fold[recipient]);
    }
  }

  free(f);
  free(fold);

  yepLibrary_Release();

  return(alpha);
}
