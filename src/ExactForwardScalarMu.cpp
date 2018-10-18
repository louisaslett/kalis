#include <Rcpp.h>
using namespace Rcpp;
#include <stdlib.h>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <functional>

#include "Cache.h"

// #include <iacaMarks.h>

void ExactForwardNoExpAVX3_scmu_cpp_raw(double *const __restrict__ alpha,
                                        double *const __restrict__ alpha_f,
                                        double *const __restrict__ alpha_f2,
                                        const int alpha_from_rec,
                                        const int alpha_t,
                                        const int t,
                                        const int from_rec,
                                        const int to_rec,
                                        const int L,
                                        const int N,
                                        const double *const __restrict__ Pi,
                                        const double mu,
                                        const double *const __restrict__ rho) {

#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)

  int_fast32_t l=alpha_t;

  double *__restrict__ f, *__restrict__ fold, *__restrict__ foldold;
  f       = (double*) malloc(sizeof(double)*(to_rec-from_rec));
  fold    = &(alpha_f[0]);
  foldold = &(alpha_f2[0]);

  // Locus zero setup
  if(l<0) {
    l = 0;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_alpha = recipient-alpha_from_rec;
      int32_t recipient_hap = (hap_locus[0][recipient/32] >> recipient%32) & 1;

      fold[recipient_alpha]    = 0.0;
      foldold[recipient_alpha] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (hap_locus[0][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu
                          + (1-H) * (1.0 - mu));

        fold[recipient_alpha] += alpha[donor + N*recipient_alpha] = theta*Pi[donor + N*recipient];
      }

      fold[recipient_alpha] = -log(fold[recipient_alpha]*rho[0]);
    }
  }

  const int_fast32_t reset_l = l;

  // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
  const double muTmp1 = 2.0 * mu - 1.0, muTmp2 = - mu + 1.0;
  const __m256d _muTmp1 = _mm256_broadcast_sd(&muTmp1), _muTmp2 = _mm256_broadcast_sd(&muTmp2);
  for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
    int_fast32_t recipient_alpha = recipient-alpha_from_rec;
    l = reset_l;

    while(l<t) {
      ++l;
      // (1.0-rho[l-1]) * fratio .... for all recipients we consider to get us going
      double fratioMulOmRho = (1.0 - rho[l-1]) * exp(fold[recipient_alpha] - foldold[recipient_alpha]);
      __m256d _fratioMulOmRho = _mm256_set1_pd(fratioMulOmRho);

      // Load this recipient's bit into all 256-bits of an AVX register
      int32_t recipient_hap = 0;
      recipient_hap -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
      __m256i _recipient_hap = _mm256_set1_epi32(recipient_hap);

      f[recipient-from_rec] = 0.0; // For accumulating scalar part in ragged end ...
      __m256d _f = _mm256_set1_pd(0.0); // ... and for vector part.

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ alphaRow = &(alpha[N*recipient_alpha]);
      const double *__restrict__ PiRow = &(Pi[N*recipient]);

      for(int_fast32_t donoroff=0; donoroff<N/(32*8); ++donoroff) {
        // Load next 256 donors and XOR with recipients
        __m256i _HA = _mm256_xor_si256(_recipient_hap, _mm256_load_si256((__m256i*) &(hap_locus[l][donoroff*8])));
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

          _f              = _mm256_add_pd(_f, _alpha1);
          _f              = _mm256_add_pd(_f, _alpha2);
          _f              = _mm256_add_pd(_f, _alpha3);
          _f              = _mm256_add_pd(_f, _alpha4);

          _mm256_storeu_pd(alphaNow1, _alpha1);
          _mm256_storeu_pd(alphaNow2, _alpha2);
          _mm256_storeu_pd(alphaNow3, _alpha3);
          _mm256_storeu_pd(alphaNow4, _alpha4);
        }
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=0; donor<N%(32*8); ++donor) {
        int32_t donor_hap = (hap_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        f[recipient-from_rec] += alphaRow[(N/(32*8))*32*8+donor] = (H * muTmp1 + muTmp2) * (PiRow[(N/(32*8))*32*8+donor] + fratioMulOmRho * alphaRow[(N/(32*8))*32*8+donor]);
      }

      // Accumulate row sum into f[recipient-from_rec]
      _f = _mm256_hadd_pd(_f, _f);
      f[recipient-from_rec] += ((double*)&_f)[0] + ((double*)&_f)[2];

      foldold[recipient_alpha] = fold[recipient_alpha];
      fold[recipient_alpha] = -(log(f[recipient-from_rec] * rho[l]) - fold[recipient_alpha]);
    }
  }

  free(f);

#endif

}

// [[Rcpp::export]]
void ExactForwardNoExpAVX3_scmu_cpp(NumericMatrix alpha,
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
                                    const double mu,
                                    NumericVector rho) {
  ExactForwardNoExpAVX3_scmu_cpp_raw(&(alpha[0]),
                                     &(alpha_f[0]),
                                     &(alpha_f2[0]),
                                     alpha_from_rec,
                                     alpha_t,
                                     t,
                                     from_rec,
                                     to_rec,
                                     L,
                                     N,
                                     &(Pi[0]),
                                     mu,
                                     &(rho[0]));
}

// [[Rcpp::export]]
void ParExactForwardNoExpAVX3_scmu_cpp(NumericMatrix alpha,
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
                                       const double mu,
                                       NumericVector rho,
                                       const int nthreads) {
  std::vector<std::thread> threads;

  if(nthreads < 2) {
    Rcout << "Only use parallel function with at least 2 threads";
  }

  // round(hap(0, nthreads) * double(to_rec-from_rec) / double(nthreads)) + from_rec;
  double spacing = double(to_rec-from_rec) / double(nthreads);

  for(int_fast32_t i=0; i<nthreads; ++i) {
    threads.push_back(std::thread(ExactForwardNoExpAVX3_scmu_cpp_raw,
                                  &(alpha[0]),
                                  &(alpha_f[0]),
                                  &(alpha_f2[0]),
                                  alpha_from_rec,
                                  alpha_t,
                                  t,
                                  round(from_rec + i*spacing),
                                  round(from_rec + (i+1)*spacing),
                                  L,
                                  N,
                                  &(Pi[0]),
                                  mu,
                                  &(rho[0])));
    // Rcout << "From: " << round(from_rec + i*spacing) << ", To: " << round(from_rec + (i+1)*spacing) << "\n";
  }

  for(auto& th : threads) {
    th.join();
  }
}
