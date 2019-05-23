#include <Rcpp.h>
using namespace Rcpp;
#include <stdlib.h>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <functional>

#include "Cache.h"

// #include <iacaMarks.h>

#ifdef EXACTBACKWARDNOEXP

#include "Stencil2.h"

void CPP_RAW_FN(EXACTBACKWARDNOEXP)(double *const __restrict__ beta,
                double *const __restrict__ beta_g,
                double *const __restrict__ beta_g2,
                const int beta_from_rec,
                const int beta_t,
                const int t,
                const int from_rec,
                const int to_rec,
                const int L,
                const int N,
                PI_TYPE_C Pi,
                const MU_TYPE_C mu,
                const double *const __restrict__ rho) {

#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)

  int_fast32_t l = beta_t;

  double *__restrict__ g, *__restrict__ gold, *__restrict__ goldold;
  g       = (double*) malloc(sizeof(double)*(to_rec-from_rec));
  gold    = &(beta_g[0]);
  goldold = &(beta_g2[0]);

  // Locus L setup
  if(l>L-1) {
    l = L-1;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
      int32_t recipient_hap = (hap_locus[l][recipient/32] >> recipient%32) & 1;

      gold[recipient_beta] = 0.0;
      goldold[recipient_beta] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (hap_locus[l][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
#if KALIS_MU == MU_SCALAR
        double theta = (H * mu
                          + (1-H) * (1.0 - mu));
#elif KALIS_MU == MU_VECTOR
        double theta = (H * mu[l]
                          + (1-H) * (1.0 - mu[l]));
#endif

        beta[donor + N*recipient_beta] = 1.0;

#if KALIS_PI == PI_SCALAR
        gold[recipient_beta] += Pi * theta;
#elif KALIS_PI == PI_MATRIX
        gold[recipient_beta] += Pi[donor + N*recipient] * theta;
#endif
      }

#if KALIS_PI == PI_SCALAR
      // Adjust due to scalar Pi for recipient/receipient locus
#if KALIS_MU == MU_SCALAR
      gold[recipient_beta] -= (1.0-mu)*Pi;
#elif KALIS_MU == MU_VECTOR
      gold[recipient_beta] -= (1.0-mu[0])*Pi;
#endif
#endif

      gold[recipient_beta] = -log(gold[recipient_beta]);
    }
  }

  const int_fast32_t reset_l = l;

#if KALIS_PI == PI_SCALAR
  const __m256d _Pi = _mm256_set1_pd(Pi);
#endif

#if KALIS_MU == MU_SCALAR
  // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
  const double muTmp1 = 2.0 * mu - 1.0, muTmp2 = - mu + 1.0;
  const __m256d _muTmp1 = _mm256_broadcast_sd(&muTmp1), _muTmp2 = _mm256_broadcast_sd(&muTmp2);
#endif

  for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
    int_fast32_t recipient_beta = recipient-beta_from_rec;
    l = reset_l;

    while(l>t) {
      --l;
      // (1.0-rho[l]) * gratio .... for all recipients we consider to get us going
      double gratioMulOmRho = (1.0 - rho[l]) * exp(gold[recipient_beta] - goldold[recipient_beta]);
      __m256d _gratioMulOmRho = _mm256_set1_pd(gratioMulOmRho);

      // Load this recipient's bit into all 256-bits of an AVX register
      int32_t recipient_hap = 0, recipient_hap_prev = 0;
      recipient_hap_prev -= (hap_locus[l+1][recipient/32] >> recipient%32) & 1;
      recipient_hap      -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
      __m256i _recipient_hap_prev = _mm256_set1_epi32(recipient_hap_prev);
      __m256i _recipient_hap      = _mm256_set1_epi32(recipient_hap);

      g[recipient-from_rec] = 0.0; // For accumulating scalar part in ragged end ...
      __m256d _g = _mm256_set1_pd(0.0); // ... and for vector part.

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ betaRow = &(beta[N*recipient_beta]);
#if KALIS_PI == PI_MATRIX
      double *__restrict__ PiRow = &(Pi[N*recipient]);
      // The diagonal should be zero ... we can fix that after the loop, but in order to avoid accumulating
      // incorrectly in g, we want to force the diagonal of Pi to zero
      PiRow[recipient] = 0.0;
#endif

#if KALIS_MU == MU_VECTOR
      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      const double muTmp1a = 2.0 * mu[l+1] - 1.0, muTmp2a = - mu[l+1] + 1.0;
      const __m256d _muTmp1a = _mm256_broadcast_sd(&muTmp1a), _muTmp2a = _mm256_broadcast_sd(&muTmp2a);
      const double muTmp1b = 2.0 * mu[l] - 1.0, muTmp2b = - mu[l] + 1.0;
      const __m256d _muTmp1b = _mm256_broadcast_sd(&muTmp1b), _muTmp2b = _mm256_broadcast_sd(&muTmp2b);
#endif

      for(int_fast32_t donoroff=0; donoroff<N/(32*8); ++donoroff) {
        // Load next 256 donors and XOR with recipients
        __m256i _HA = _mm256_xor_si256(_recipient_hap_prev, _mm256_load_si256((__m256i*) &(hap_locus[l+1][donoroff*8])));
        uint32_t *HA = (uint32_t*) &_HA;
        __m256i _HB = _mm256_xor_si256(_recipient_hap,      _mm256_load_si256((__m256i*) &(hap_locus[l][donoroff*8])));
        uint32_t *HB = (uint32_t*) &_HB;

        const uint32_t mask = 16843009;
        for(int_fast32_t donor=0; donor<((32*8)/4)/4; ++donor) {
          // IACA_START
          double *betaNow1 = betaRow + donoroff*32*8 + donor*4*4;
          double *betaNow2 = betaRow + donoroff*32*8 + donor*4*4 + 4;
          double *betaNow3 = betaRow + donoroff*32*8 + donor*4*4 + 8;
          double *betaNow4 = betaRow + donoroff*32*8 + donor*4*4 + 12;

          __m256d _Rho = _mm256_set1_pd(rho[l]);

          __m256d _beta1  = _mm256_loadu_pd(betaNow1);
          __m256d _beta2  = _mm256_loadu_pd(betaNow2);
          __m256d _beta3  = _mm256_loadu_pd(betaNow3);
          __m256d _beta4  = _mm256_loadu_pd(betaNow4);

          __m256d _theta1 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4)/8]) >> (((donor*4)%8)*4), mask))));
          __m256d _theta2 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+1)/8]) >> (((donor*4+1)%8)*4), mask))));
          __m256d _theta3 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+2)/8]) >> (((donor*4+2)%8)*4), mask))));
          __m256d _theta4 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+3)/8]) >> (((donor*4+3)%8)*4), mask))));

#if KALIS_MU == MU_SCALAR
          _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#elif KALIS_MU == MU_VECTOR
          _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
          _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
          _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
          _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
#endif

          _beta1          = _mm256_mul_pd(_beta1, _theta1); // (theta*beta)
          _beta2          = _mm256_mul_pd(_beta2, _theta2); // (theta*beta)
          _beta3          = _mm256_mul_pd(_beta3, _theta3); // (theta*beta)
          _beta4          = _mm256_mul_pd(_beta4, _theta4); // (theta*beta)

          _beta1          = _mm256_fmadd_pd(_beta1, _gratioMulOmRho, _Rho); // (rho + {theta*beta} * {(1-rho)gratio})
          _beta2          = _mm256_fmadd_pd(_beta2, _gratioMulOmRho, _Rho); // (rho + {theta*beta} * {(1-rho)gratio})
          _beta3          = _mm256_fmadd_pd(_beta3, _gratioMulOmRho, _Rho); // (rho + {theta*beta} * {(1-rho)gratio})
          _beta4          = _mm256_fmadd_pd(_beta4, _gratioMulOmRho, _Rho); // (rho + {theta*beta} * {(1-rho)gratio})

#if KALIS_PI == PI_MATRIX
          __m256d _pi1    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4);
          __m256d _pi2    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 4);
          __m256d _pi3    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 8);
          __m256d _pi4    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 12);
#endif

          _theta1         = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HB[(donor*4)/8]) >> (((donor*4)%8)*4), mask))));
          _theta2         = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HB[(donor*4+1)/8]) >> (((donor*4+1)%8)*4), mask))));
          _theta3         = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HB[(donor*4+2)/8]) >> (((donor*4+2)%8)*4), mask))));
          _theta4         = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HB[(donor*4+3)/8]) >> (((donor*4+3)%8)*4), mask))));

#if KALIS_MU == MU_SCALAR
          _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#elif KALIS_MU == MU_VECTOR
          _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
          _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
          _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
          _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
#endif

          _theta1         = _mm256_mul_pd(_beta1, _theta1); // (theta*beta)
          _theta2         = _mm256_mul_pd(_beta2, _theta2); // (theta*beta)
          _theta3         = _mm256_mul_pd(_beta3, _theta3); // (theta*beta)
          _theta4         = _mm256_mul_pd(_beta4, _theta4); // (theta*beta)

#if KALIS_PI == PI_MATRIX
          _g              = _mm256_fmadd_pd(_pi1, _theta1, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_pi2, _theta2, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_pi3, _theta3, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_pi4, _theta4, _g); // g += Pi * {theta*beta}
#elif KALIS_PI == PI_SCALAR
          _g              = _mm256_fmadd_pd(_Pi, _theta1, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_Pi, _theta2, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_Pi, _theta3, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_Pi, _theta4, _g); // g += Pi * {theta*beta}
#endif

          _mm256_storeu_pd(betaNow1, _beta1);
          _mm256_storeu_pd(betaNow2, _beta2);
          _mm256_storeu_pd(betaNow3, _beta3);
          _mm256_storeu_pd(betaNow4, _beta4);
        }
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=0; donor<N%(32*8); ++donor) {
        int32_t donor_hap = (hap_locus[l+1][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
        int32_t H = (recipient_hap_prev ^ donor_hap) & 1;
#if KALIS_MU == MU_SCALAR
        betaRow[(N/(32*8))*32*8+donor] = rho[l] + (H * muTmp1 + muTmp2) * betaRow[(N/(32*8))*32*8+donor] * gratioMulOmRho;
#elif KALIS_MU == MU_VECTOR
        betaRow[(N/(32*8))*32*8+donor] = rho[l] + (H * muTmp1a + muTmp2a) * betaRow[(N/(32*8))*32*8+donor] * gratioMulOmRho;
#endif

        donor_hap = (hap_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
        H = (recipient_hap ^ donor_hap) & 1;
#if KALIS_PI == PI_MATRIX && KALIS_MU == MU_VECTOR
        g[recipient-from_rec] += PiRow[(N/(32*8))*32*8+donor] * (H * muTmp1b + muTmp2b) * betaRow[(N/(32*8))*32*8+donor];
#elif KALIS_PI == PI_MATRIX && KALIS_MU == MU_SCALAR
        g[recipient-from_rec] += PiRow[(N/(32*8))*32*8+donor] * (H * muTmp1 + muTmp2) * betaRow[(N/(32*8))*32*8+donor];
#elif KALIS_PI == PI_SCALAR && KALIS_MU == MU_VECTOR
        g[recipient-from_rec] += Pi * (H * muTmp1b + muTmp2b) * betaRow[(N/(32*8))*32*8+donor];
#elif KALIS_PI == PI_SCALAR && KALIS_MU == MU_SCALAR
        g[recipient-from_rec] += Pi * (H * muTmp1 + muTmp2) * betaRow[(N/(32*8))*32*8+donor];
#endif
      }
#if KALIS_PI == PI_SCALAR
      // Adjustments for scalar Pi
#if KALIS_MU == MU_VECTOR
      g[recipient-from_rec] -= Pi * muTmp2b * betaRow[recipient];
#elif KALIS_MU == MU_SCALAR
      g[recipient-from_rec] -= Pi * muTmp2 * betaRow[recipient];
#endif
#elif KALIS_PI == PI_MATRIX
      // Fix the diagonal of beta ... note that because we set the diagonal of Pi
      // to zero we don't need to adjust g
#endif
      betaRow[recipient] = 0.0;

      // Accumulate row sum into g[recipient-from_rec]
      _g = _mm256_hadd_pd(_g, _g);
      g[recipient-from_rec] += ((double*)&_g)[0] + ((double*)&_g)[2];

      goldold[recipient_beta] = gold[recipient_beta];
      gold[recipient_beta] = -(log(g[recipient-from_rec]) - gold[recipient_beta]);
    }
  }

  free(g);

#endif

}



void CPP_FN(EXACTBACKWARDNOEXP)(NumericMatrix beta,
            NumericVector beta_g,
            NumericVector beta_g2,
            const int beta_from_rec,
            const int beta_t,
            const int t,
            const int from_rec,
            const int to_rec,
            const int L,
            const int N,
            PI_TYPE_CPP Pi,
            MU_TYPE_CPP mu,
            NumericVector rho) {
  CPP_RAW_FN(EXACTBACKWARDNOEXP)(&(beta[0]),
             &(beta_g[0]),
             &(beta_g2[0]),
             beta_from_rec,
             beta_t,
             t,
             from_rec,
             to_rec,
             L,
             N,
             PI_ARG_CPP,
             MU_ARG_CPP,
             &(rho[0]));
}



void PAR_CPP_FN(EXACTBACKWARDNOEXP)(NumericMatrix beta,
                NumericVector beta_g,
                NumericVector beta_g2,
                const int beta_from_rec,
                const int beta_t,
                const int t,
                const int from_rec,
                const int to_rec,
                const int L,
                const int N,
                PI_TYPE_CPP Pi,
                MU_TYPE_CPP mu,
                NumericVector rho,
                const int nthreads) {
  std::vector<std::thread> threads;

  if(nthreads < 2) {
    Rcout << "Only use parallel function with at least 2 threads";
  }

  // round(hap(0, nthreads) * double(to_rec-from_rec) / double(nthreads)) + from_rec;
  double spacing = double(to_rec-from_rec) / double(nthreads);

  for(int_fast32_t i=0; i<nthreads; ++i) {
    threads.push_back(std::thread(CPP_RAW_FN(EXACTBACKWARDNOEXP),
                                  &(beta[0]),
                                  &(beta_g[0]),
                                  &(beta_g2[0]),
                                  beta_from_rec,
                                  beta_t,
                                  t,
                                  round(from_rec + i*spacing),
                                  round(from_rec + (i+1)*spacing),
                                  L,
                                  N,
                                  PI_ARG_CPP,
                                  MU_ARG_CPP,
                                  &(rho[0])));
    // Rcout << "From: " << round(from_rec + i*spacing) << ", To: " << round(from_rec + (i+1)*spacing) << "\n";
  }

  for(auto& th : threads) {
    th.join();
  }
}

#endif
