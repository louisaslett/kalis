#include <Rcpp.h>
using namespace Rcpp;

#include <stdlib.h>
#include <yepLibrary.h>
#include <yepCore.h>
#include <yepMath.h>

#include <immintrin.h>

#include <iacaMarks.h>

#include "Cache.h"

// [[Rcpp::export]]
void ExactForwardYepppExpC_cpp(NumericMatrix alpha, NumericVector alpha_f, int alpha_t, int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
  int_fast32_t l=alpha_t;

  yepLibrary_Init();

  double *f, *fold;
  f    = (double*) malloc(sizeof(double)*N);
  fold = &(alpha_f[0]);

  // Locus zero setup
  if(l<0) {
    l = 0;
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

  yepLibrary_Release();
}

// [[Rcpp::export]]
void ExactForwardYepppExpAVX_cpp(NumericMatrix alpha, NumericVector alpha_f, int alpha_t, int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
  int_fast32_t l=alpha_t;

  yepLibrary_Init();

  double *f, *fold;
  f    = (double*) malloc(sizeof(double)*N);
  fold = &(alpha_f[0]);

  // Locus zero setup
  if(l<0) {
    l = 0;
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
  }

  while(l<t) {
    ++l;
    for(int_fast32_t recipient=0; recipient<N; ++recipient) {
      int32_t recipient_hap = 0;
      recipient_hap -= (seq_locus[l][recipient/32] >> recipient%32) & 1;
      __m256i _recipient_hap = _mm256_set1_epi32(recipient_hap);

      f[recipient] = 0.0;

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *alphaRow = &(alpha[N*recipient]);
      double YepppTmp[N] __attribute__ ((aligned (32)));
      yepCore_Add_V64fS64f_V64f(alphaRow, fold[recipient], YepppTmp, N);
      yepMath_Exp_V64f_V64f(YepppTmp, alphaRow, N);

      // (1-rho) * alphaRow[donor]
      yepCore_Multiply_IV64fS64f_IV64f(alphaRow, (1.0-rho[l-1]), N);
      // theta[donor] * ((1-rho) * alphaRow[donor])
      double muTmp1 = 2.0 * mu[l] - 1.0, muTmp2 = - mu[l] + 1.0;
      __m256d _muTmp1 = _mm256_broadcast_sd(&muTmp1), _muTmp2 = _mm256_broadcast_sd(&muTmp2);
      for(int_fast32_t donoroff=0; donoroff<N/(32*8); ++donoroff) {
        // Tried pipelining two packed doubles but it's slower (register pressure?)
        __m256i _donor_hapA = _mm256_load_si256((__m256i*) &(seq_locus[l][donoroff*8]));

        //int32_t H = (recipient_hap ^ seq_locus[l][donoroff]); // _mm256_xor_si256
        __m256i _HA = _mm256_xor_si256(_recipient_hap, _donor_hapA);
        __m128i ones = _mm_set1_epi32(1);
        double tmp[8] __attribute__ ((aligned (32)));
        for(char donor=0; donor<32; ++donor) {
        // donor_hap = (seq_locus[l][donoroff] >> donor) & 1;
        // H = (recipient_hap ^ donor_hap) & 1;
        // YepppTmp[donoroff*8+donor] = (H * muTmp
        //                      + (1-H) * (1.0 - muTmp));
        // Extract the low and high 128-bits
        __m128i _HAlo, _HAhi;

        //_mm256_storeu2_m128i(&_HAhi, &_HAlo, _HA); // Clang only
        _HAlo = _mm256_castsi256_si128(_HA);
        _HAhi = _mm256_extractf128_si256(_HA, 1);

        _HAlo = _mm_and_si128(_HAlo, ones);
        _HAhi = _mm_and_si128(_HAhi, ones);

        // Convert to packed doubles and do the theta computations
        _mm256_store_pd(tmp,   _mm256_fmadd_pd(_mm256_cvtepi32_pd(_HAlo), _muTmp1, _muTmp2)); // YepppTmp[donoroff*32+donor] = (H&1) * (2.0 * mu[l] - 1.0) - mu[l] + 1.0;

        // Scatter to memory -- sadly scatter intrinsics are only in AVX-512
        YepppTmp[donoroff*32*8       +donor] = tmp[0];
        YepppTmp[donoroff*32*8 +32   +donor] = tmp[1];
        __m256d __HAhi = _mm256_cvtepi32_pd(_HAhi);
        YepppTmp[donoroff*32*8 +32*2 +donor] = tmp[2];
        YepppTmp[donoroff*32*8 +32*3 +donor] = tmp[3];
        _mm256_store_pd(tmp+4, _mm256_fmadd_pd(__HAhi, _muTmp1, _muTmp2)); // YepppTmp[donoroff*32+donor] = (H&1) * (2.0 * mu[l] - 1.0) - mu[l] + 1.0;
        YepppTmp[donoroff*32*8 +32*4 +donor] = tmp[4];
        YepppTmp[donoroff*32*8 +32*5 +donor] = tmp[5];
        YepppTmp[donoroff*32*8 +32*6 +donor] = tmp[6];
        YepppTmp[donoroff*32*8 +32*7 +donor] = tmp[7];

        // Rotate to next individual in all slots
        _HA = _mm256_srli_epi32(_HA, 1); // H >>= 1;
        }
      }
        // Tidy up any ragged end past a multiple of 256 ...
        for(int32_t donor=0; donor<N%(32*8); ++donor) {
        int32_t donor_hap = (seq_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        YepppTmp[(N/(32*8))*32*8+donor] = H * muTmp1 + muTmp2;
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

        yepLibrary_Release();
}

      // [[Rcpp::export]]
      void ExactForwardYepppExpAVX2_cpp(NumericMatrix alpha, NumericVector alpha_f, int alpha_t, int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
      int_fast32_t l=alpha_t;

      yepLibrary_Init();

      double *f, *fold;
      f    = (double*) malloc(sizeof(double)*N);
      fold = &(alpha_f[0]);

      // Locus zero setup
      if(l<0) {
      l = 0;
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
      }

      while(l<t) {
      ++l;
      for(int_fast32_t recipient=0; recipient<N; ++recipient) {
      int32_t recipient_hap = 0;
      recipient_hap -= (seq_locus[l][recipient/32] >> recipient%32) & 1;
      __m256i _recipient_hap = _mm256_set1_epi32(recipient_hap);

      f[recipient] = 0.0;

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *alphaRow = &(alpha[N*recipient]);
      double YepppTmp[N] __attribute__ ((aligned (32)));
      yepCore_Add_V64fS64f_V64f(alphaRow, fold[recipient], YepppTmp, N);
      yepMath_Exp_V64f_V64f(YepppTmp, alphaRow, N);

      // (1-rho) * alphaRow[donor]
      yepCore_Multiply_IV64fS64f_IV64f(alphaRow, (1.0-rho[l-1]), N);
      // theta[donor] * ((1-rho) * alphaRow[donor])
      double muTmp1 = 2.0 * mu[l] - 1.0, muTmp2 = - mu[l] + 1.0;
      __m256d _muTmp1 = _mm256_broadcast_sd(&muTmp1), _muTmp2 = _mm256_broadcast_sd(&muTmp2);
      for(int_fast32_t donoroff=0; donoroff<N/(32*8); ++donoroff) {
      // Tried pipelining two packed doubles but it's slower (register pressure?)
      __m256i _donor_hapA = _mm256_load_si256((__m256i*) &(seq_locus[l][donoroff*8]));

      //int32_t H = (recipient_hap ^ seq_locus[l][donoroff]); // _mm256_xor_si256
      __m256i _HA = _mm256_xor_si256(_recipient_hap, _donor_hapA);
      __m128i ones = _mm_set1_epi32(1);
      double tmp[8] __attribute__ ((aligned (32)));
      for(char donor=0; donor<32; ++donor) {
        // donor_hap = (seq_locus[l][donoroff] >> donor) & 1;
        // H = (recipient_hap ^ donor_hap) & 1;
        // YepppTmp[donoroff*8+donor] = (H * muTmp
                                         //                      + (1-H) * (1.0 - muTmp));
        // Extract the low and high 128-bits
        __m128i _HAlo, _HAhi;

        //_mm256_storeu2_m128i(&_HAhi, &_HAlo, _HA); // Clang only
        _HAlo = _mm256_castsi256_si128(_HA);
        _HAhi = _mm256_extractf128_si256(_HA, 1);

        _HAlo = _mm_and_si128(_HAlo, ones);
        _HAhi = _mm_and_si128(_HAhi, ones);

        // Convert to packed doubles and do the theta computations
        _mm256_store_pd(tmp,   _mm256_fmadd_pd(_mm256_cvtepi32_pd(_HAlo), _muTmp1, _muTmp2)); // YepppTmp[donoroff*32+donor] = (H&1) * (2.0 * mu[l] - 1.0) - mu[l] + 1.0;
        _mm256_store_pd(tmp+4, _mm256_fmadd_pd(_mm256_cvtepi32_pd(_HAhi), _muTmp1, _muTmp2)); // YepppTmp[donoroff*32+donor] = (H&1) * (2.0 * mu[l] - 1.0) - mu[l] + 1.0;

        // Scatter to memory -- sadly scatter intrinsics are only in AVX-512
        YepppTmp[donoroff*32*8       +donor] = tmp[0];
        YepppTmp[donoroff*32*8 +32   +donor] = tmp[1];
        YepppTmp[donoroff*32*8 +32*2 +donor] = tmp[2];
        YepppTmp[donoroff*32*8 +32*3 +donor] = tmp[3];
        YepppTmp[donoroff*32*8 +32*4 +donor] = tmp[4];
        YepppTmp[donoroff*32*8 +32*5 +donor] = tmp[5];
        YepppTmp[donoroff*32*8 +32*6 +donor] = tmp[6];
        YepppTmp[donoroff*32*8 +32*7 +donor] = tmp[7];

        // Rotate to next individual in all slots
        _HA = _mm256_srli_epi32(_HA, 1); // H >>= 1;
      }
      }
    // Tidy up any ragged end past a multiple of 256 ...
    for(int32_t donor=0; donor<N%(32*8); ++donor) {
      int32_t donor_hap = (seq_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
      int32_t H = (recipient_hap ^ donor_hap) & 1;
      YepppTmp[(N/(32*8))*32*8+donor] = H * muTmp1 + muTmp2;
    }
    yepCore_Multiply_IV64fV64f_IV64f(alphaRow, YepppTmp, N);

    // TODO: benchmark against a fused multiply and add for these two together
    // theta[donor] * Pi(donor, recipient)
    //yepCore_Multiply_IV64fV64f_IV64f(YepppTmp, &(Pi(0, recipient)), N);
    // Tot up alphaRow
    //yepCore_Add_IV64fV64f_IV64f(alphaRow, YepppTmp, N);
    double *PiRow = &(Pi(0, recipient));
    for(int_fast32_t donor=0; donor<N/(4*4); ++donor) {
      // _mm256_storeu_pd(alphaRow + donor*4*2,     _mm256_fmadd_pd(_mm256_loadu_pd(YepppTmp + donor*4*2),
                                                                     //                                                           _mm256_loadu_pd(PiRow + donor*4*2),
                                                                     //                                                           _mm256_loadu_pd(alphaRow + donor*4*2)));
      // _mm256_storeu_pd(alphaRow + donor*4*2 + 4, _mm256_fmadd_pd(_mm256_loadu_pd(YepppTmp + donor*4*2 + 4),
                                                                     //                                                           _mm256_loadu_pd(PiRow + donor*4*2 + 4),
                                                                     //                                                           _mm256_loadu_pd(alphaRow + donor*4*2 + 4)));
      __m256d A1 = _mm256_loadu_pd(YepppTmp + donor*4*4);
      __m256d B1 = _mm256_loadu_pd(PiRow + donor*4*4);
      __m256d A2 = _mm256_loadu_pd(YepppTmp + donor*4*4 + 4);
      __m256d B2 = _mm256_loadu_pd(PiRow + donor*4*4 + 4);
      __m256d A3 = _mm256_loadu_pd(YepppTmp + donor*4*4 + 8);
      __m256d B3 = _mm256_loadu_pd(PiRow + donor*4*4 + 8);
      __m256d A4 = _mm256_loadu_pd(YepppTmp + donor*4*4 + 12);
      __m256d B4 = _mm256_loadu_pd(PiRow + donor*4*4 + 12);
      _mm256_storeu_pd(alphaRow + donor*4*4, _mm256_fmadd_pd(A1, B1, _mm256_loadu_pd(alphaRow + donor*4*4)));
      _mm256_storeu_pd(alphaRow + donor*4*4 + 4, _mm256_fmadd_pd(A2, B2, _mm256_loadu_pd(alphaRow + donor*4*4 + 4)));
      _mm256_storeu_pd(alphaRow + donor*4*4 + 8, _mm256_fmadd_pd(A3, B3, _mm256_loadu_pd(alphaRow + donor*4*4 + 8)));
      _mm256_storeu_pd(alphaRow + donor*4*4 + 12, _mm256_fmadd_pd(A4, B4, _mm256_loadu_pd(alphaRow + donor*4*4 + 12)));
    }
    for(int_fast32_t donor=0; donor<N%(4*4); ++donor) {
      alphaRow[(N/(4*4))*4*4 + donor] += YepppTmp[(N/(4*4))*4*4 + donor]*Pi((N/(4*4))*4*4 + donor, recipient);
    }

    yepCore_Sum_V64f_S64f(alphaRow, f+recipient, N);
    // f[recipient] += alphaRow[donor] = (theta * PiRow[donor]
                                          //                                      + theta * (1.0-rhoTmp[0]) * alphaRow[donor]);

    yepMath_Log_V64f_V64f(alphaRow, YepppTmp, N);
    yepCore_Subtract_V64fS64f_V64f(YepppTmp, fold[recipient], alphaRow, N);

    fold[recipient] = -(log(f[recipient] * rho[l]) - fold[recipient]);
      }
      }

free(f);

yepLibrary_Release();
      }

// [[Rcpp::export]]
void ExactForwardNoExpAVX_cpp(NumericMatrix alpha, NumericVector alpha_f, NumericVector alpha_f2, int alpha_t, int t, int from_rec, int to_rec, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
  int_fast32_t l=alpha_t;

  yepLibrary_Init();

  double *__restrict__ f, *__restrict__ fold, *__restrict__ foldold;
  f    = (double*) malloc(sizeof(double)*N);
  fold = &(alpha_f[0]);
  foldold = &(alpha_f2[0]);

  // Locus zero setup
  if(l<0) {
    l = 0;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int32_t recipient_hap = (seq_locus[0][recipient/32] >> recipient%32) & 1;
      fold[recipient]    = 0.0;
      foldold[recipient] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (seq_locus[0][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu[0]
                        + (1-H) * (1.0 - mu[0]));

        fold[recipient] += alpha(donor, recipient) = theta*Pi(donor, recipient);

        // alpha(donor, recipient) = log(alpha(donor, recipient));
      }

      fold[recipient] = -log(fold[recipient]*rho[0]);
    }
  }

  while(l<t) {
    ++l;

    // (1.0-rho[l-1]) * fratio .... for all recipients we consider to get us going
    double fratioMulOmRho[to_rec-from_rec] __attribute__ ((aligned (64)));
    yepCore_Subtract_V64fIV64f_IV64f(fold+from_rec, foldold+from_rec, to_rec-from_rec);
    yepMath_Exp_V64f_V64f(foldold+from_rec, fratioMulOmRho, to_rec-from_rec);
    yepCore_Multiply_IV64fS64f_IV64f(fratioMulOmRho, 1.0-rho[l-1], to_rec-from_rec);

    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      // Load this recipient's bit into all 256-bits of an AVX register
      int32_t recipient_hap = 0;
      recipient_hap -= (seq_locus[l][recipient/32] >> recipient%32) & 1;
      __m256i _recipient_hap = _mm256_set1_epi32(recipient_hap);

      f[recipient] = 0.0;

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ alphaRow = &(alpha[N*recipient]);
      double *__restrict__ PiRow = &(Pi[N*recipient]);
      // Load (1.0-rho[l-1]) * fratio into AVX register
      __m256d _fratioMulOmRho = _mm256_broadcast_sd(fratioMulOmRho+recipient-from_rec);

      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      double muTmp1 = 2.0 * mu[l] - 1.0, muTmp2 = - mu[l] + 1.0;
      __m256d _muTmp1 = _mm256_broadcast_sd(&muTmp1), _muTmp2 = _mm256_broadcast_sd(&muTmp2);
      for(int_fast32_t donoroff=0; donoroff<N/(32*8); ++donoroff) {
      // Load next 256 donors and XOR with recipients
      __m256i _HA = _mm256_xor_si256(_recipient_hap, _mm256_load_si256((__m256i*) &(seq_locus[l][donoroff*8])));
      uint32_t *HA = (uint32_t*) &_HA;

      const uint32_t mask = 16843009;
      for(int_fast32_t donor=0; donor<(32*8)/4; ++donor) {
      // IACA_START
      __m256d _alpha1 = _mm256_loadu_pd(alphaRow + donoroff*32*8 + donor*4);
      __m256d _pi1    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4);
      _alpha1         = _mm256_fmadd_pd(_alpha1, _fratioMulOmRho, _pi1); // (Pi + {(1-rho)*fratio} * alpha)
      __m256d _theta1 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[donor/8]) >> ((donor%8)*4), mask))));
      _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
      _alpha1         = _mm256_mul_pd(_theta1, _alpha1);
      _mm256_storeu_pd(alphaRow + donoroff*32*8 + donor*4, _alpha1);
      }
      // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=0; donor<N%(32*8); ++donor) {
      int32_t donor_hap = (seq_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
      int32_t H = (recipient_hap ^ donor_hap) & 1;
      alphaRow[(N/(32*8))*32*8+donor] = (H * muTmp1 + muTmp2) * (PiRow[(N/(32*8))*32*8+donor] + fratioMulOmRho[recipient-from_rec] * alphaRow[(N/(32*8))*32*8+donor]);
      }

      // Accumulate row sum into f[recipient]
      yepCore_Sum_V64f_S64f(alphaRow, f+recipient, N);

      foldold[recipient] = fold[recipient];
      fold[recipient] = -(log(f[recipient] * rho[l]) - fold[recipient]);
    }
  }

      free(f);

      yepLibrary_Release();
}

// [[Rcpp::export]]
void ExactForwardNoExpAVX2_cpp(NumericMatrix alpha,
                               NumericVector alpha_f,
                               NumericVector alpha_f2,
                               const int alpha_from_rec,
                               const int alpha_t,
                               const int t,
                               const int from_rec,
                               const int to_rec,
                               const int L,
                               const int N,
                               NumericMatrix Pi,
                               NumericVector mu,
                               NumericVector rho) {
  int_fast32_t l=alpha_t;

  yepLibrary_Init();

  double *__restrict__ f, *__restrict__ fold, *__restrict__ foldold;
  f    = (double*) malloc(sizeof(double)*(to_rec-from_rec));
  fold = &(alpha_f[0]);
  foldold = &(alpha_f2[0]);

  // Locus zero setup
  if(l<0) {
    l = 0;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_alpha = recipient-alpha_from_rec;
      int32_t recipient_hap = (seq_locus[0][recipient/32] >> recipient%32) & 1;
      fold[recipient_alpha]    = 0.0;
      foldold[recipient_alpha] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (seq_locus[0][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu[0]
                          + (1-H) * (1.0 - mu[0]));

        fold[recipient_alpha] += alpha(donor, recipient_alpha) = theta*Pi(donor, recipient);
      }

      fold[recipient_alpha] = -log(fold[recipient_alpha]*rho[0]);
    }
  }

  while(l<t) {
    ++l;

    // (1.0-rho[l-1]) * fratio .... for all recipients we consider to get us going
    double fratioMulOmRho[to_rec-from_rec] __attribute__ ((aligned (64)));
    yepCore_Subtract_V64fIV64f_IV64f(fold+from_rec-alpha_from_rec, foldold+from_rec-alpha_from_rec, to_rec-from_rec);
    yepMath_Exp_V64f_V64f(foldold+from_rec-alpha_from_rec, fratioMulOmRho, to_rec-from_rec);
    yepCore_Multiply_IV64fS64f_IV64f(fratioMulOmRho, 1.0-rho[l-1], to_rec-from_rec);

    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_alpha = recipient-alpha_from_rec;

      // Load this recipient's bit into all 256-bits of an AVX register
      int32_t recipient_hap = 0;
      recipient_hap -= (seq_locus[l][recipient/32] >> recipient%32) & 1;
      __m256i _recipient_hap = _mm256_set1_epi32(recipient_hap);

      f[recipient-from_rec] = 0.0;

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ alphaRow = &(alpha[N*recipient_alpha]);
      double *__restrict__ PiRow = &(Pi[N*recipient]);
      // Load (1.0-rho[l-1]) * fratio into AVX register
      __m256d _fratioMulOmRho = _mm256_broadcast_sd(fratioMulOmRho+recipient-from_rec);

      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      double muTmp1 = 2.0 * mu[l] - 1.0, muTmp2 = - mu[l] + 1.0;
      __m256d _muTmp1 = _mm256_broadcast_sd(&muTmp1), _muTmp2 = _mm256_broadcast_sd(&muTmp2);
      for(int_fast32_t donoroff=0; donoroff<N/(32*8); ++donoroff) {
        // Load next 256 donors and XOR with recipients
        __m256i _HA = _mm256_xor_si256(_recipient_hap, _mm256_load_si256((__m256i*) &(seq_locus[l][donoroff*8])));
        uint32_t *HA = (uint32_t*) &_HA;

        const uint32_t mask = 16843009;
        for(int_fast32_t donor=0; donor<((32*8)/4)/4; ++donor) {
          // IACA_START
          double *alphaNow1 = alphaRow + donoroff*32*8 + donor*4*4;
          double *alphaNow2 = alphaRow + donoroff*32*8 + donor*4*4 + 4;
          double *alphaNow3 = alphaRow + donoroff*32*8 + donor*4*4 + 8;
          double *alphaNow4 = alphaRow + donoroff*32*8 + donor*4*4 + 12;

          __m256d _alpha1 = _mm256_loadu_pd(alphaNow1);
          __m256d _alpha2 = _mm256_loadu_pd(alphaNow2);
          __m256d _alpha3 = _mm256_loadu_pd(alphaNow3);
          __m256d _alpha4 = _mm256_loadu_pd(alphaNow4);

          __m256d _pi1    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4);
          __m256d _pi2    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 4);
          __m256d _pi3    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 8);
          __m256d _pi4    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 12);

          _alpha1         = _mm256_fmadd_pd(_alpha1, _fratioMulOmRho, _pi1); // (Pi + {(1-rho)*fratio} * alpha)
          _alpha2         = _mm256_fmadd_pd(_alpha2, _fratioMulOmRho, _pi2); // (Pi + {(1-rho)*fratio} * alpha)
          _alpha3         = _mm256_fmadd_pd(_alpha3, _fratioMulOmRho, _pi3); // (Pi + {(1-rho)*fratio} * alpha)
          _alpha4         = _mm256_fmadd_pd(_alpha4, _fratioMulOmRho, _pi4); // (Pi + {(1-rho)*fratio} * alpha)

          __m256d _theta1 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4)/8]) >> (((donor*4)%8)*4), mask))));
          __m256d _theta2 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+1)/8]) >> (((donor*4+1)%8)*4), mask))));
          __m256d _theta3 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+2)/8]) >> (((donor*4+2)%8)*4), mask))));
          __m256d _theta4 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+3)/8]) >> (((donor*4+3)%8)*4), mask))));

          _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1

          _alpha1         = _mm256_mul_pd(_theta1, _alpha1);
          _alpha2         = _mm256_mul_pd(_theta2, _alpha2);
          _alpha3         = _mm256_mul_pd(_theta3, _alpha3);
          _alpha4         = _mm256_mul_pd(_theta4, _alpha4);

          _mm256_storeu_pd(alphaNow1, _alpha1);
          _mm256_storeu_pd(alphaNow2, _alpha2);
          _mm256_storeu_pd(alphaNow3, _alpha3);
          _mm256_storeu_pd(alphaNow4, _alpha4);
        }
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=0; donor<N%(32*8); ++donor) {
        int32_t donor_hap = (seq_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        alphaRow[(N/(32*8))*32*8+donor] = (H * muTmp1 + muTmp2) * (PiRow[(N/(32*8))*32*8+donor] + fratioMulOmRho[recipient-from_rec] * alphaRow[(N/(32*8))*32*8+donor]);
      }

      // Accumulate row sum into f[recipient-from_rec]
      yepCore_Sum_V64f_S64f(alphaRow, f+recipient-from_rec, N);

      foldold[recipient_alpha] = fold[recipient_alpha];
      fold[recipient_alpha] = -(log(f[recipient-from_rec] * rho[l]) - fold[recipient_alpha]);
    }
  }

  free(f);

  yepLibrary_Release();
}
