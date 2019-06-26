#include <Rcpp.h>
using namespace Rcpp;
#include <stdlib.h>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <functional>

#include "Cache.h"

// #include <iacaMarks.h>

#ifdef EXACTFORWARDNOEXP

#include "Stencil2.h"

void CPP_RAW_FN(EXACTFORWARDNOEXP)(double *const __restrict__ alpha,
                double *const __restrict__ alpha_f,
                double *const __restrict__ alpha_f2,
                const int alpha_from_rec,
                const int alpha_t,
                const int t,
                const int from_rec,
                const int to_rec,
                const int L,
                const int N,
                const PI_TYPE_C Pi,
                const MU_TYPE_C mu,
                const double *const __restrict__ rho) {

#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)

  int_fast32_t l = alpha_t;

  double f;
  double *__restrict__ fold;
  fold    = &(alpha_f[0]);

  // Locus zero setup
  if(l<0) {
    l = 0;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_alpha = recipient-alpha_from_rec;
      int32_t recipient_hap = (hap_locus[0][recipient/32] >> recipient%32) & 1;
      recipient_hap = 1-recipient_hap; // So H below will now be 1-H

      fold[recipient_alpha]    = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (hap_locus[0][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
#if KALIS_MU == MU_SCALAR
        double theta = (H * (1.0 - 2.0*mu) + mu);
#elif KALIS_MU == MU_VECTOR
        double theta = (H * (1.0 - 2.0*mu[0]) + mu[0]);
#endif

#if KALIS_PI == PI_SCALAR
        fold[recipient_alpha] += alpha[donor + N*recipient_alpha] = theta*Pi;
#elif KALIS_PI == PI_MATRIX
        fold[recipient_alpha] += alpha[donor + N*recipient_alpha] = theta*Pi[donor + N*recipient];
#endif
      }

#if KALIS_PI == PI_SCALAR
      // Adjust due to scalar Pi for recipient/recipient locus
      // We can do this at the first locus but not later due to
      // numerical stability problems if (1-mu)*Pi is large
      // compared to all other contributions
      alpha[recipient + N*recipient_alpha] = 0.0;
#if KALIS_MU == MU_SCALAR
      fold[recipient_alpha] -= (1.0-mu)*Pi;
#elif KALIS_MU == MU_VECTOR
      fold[recipient_alpha] -= (1.0-mu[0])*Pi;
#endif
#endif
    }
  }

  const int_fast32_t reset_l = l;


#if KALIS_MU == MU_SCALAR && !defined(KALIS_1STEP)
  // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
  // ... now update to (1-H) * (1-2*mu) + mu
  const double muTmp1 = 1.0 - 2.0 * mu, muTmp2 = mu;
  const __m256d _muTmp1 = _mm256_broadcast_sd(&muTmp1), _muTmp2 = _mm256_broadcast_sd(&muTmp2);
#endif

  for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
    int_fast32_t recipient_alpha = recipient-alpha_from_rec;
#if KALIS_PI == PI_SCALAR
    const __m256d _recipient = _mm256_set1_pd(recipient);
    const __m256d _donorinc = _mm256_set1_pd(16.0);
    __m256d _donornum1 = _mm256_setr_pd(0.0, 1.0, 2.0, 3.0);
    __m256d _donornum2 = _mm256_setr_pd(4.0, 5.0, 6.0, 7.0);
    __m256d _donornum3 = _mm256_setr_pd(8.0, 9.0, 10.0, 11.0);
    __m256d _donornum4 = _mm256_setr_pd(12.0, 13.0, 14.0, 15.0);
#endif
    l = reset_l;

    while(l<t) {
      ++l;
      // (1.0-rho[l-1]) * fratio .... for all recipients we consider to get us going
      double omRhoDivF = (1.0 - rho[l-1]) / fold[recipient_alpha];
      __m256d _omRhoDivF = _mm256_set1_pd(omRhoDivF);

#if !defined(KALIS_1STEP)
      // Load this recipient's bit into all 256-bits of an AVX register
      int32_t recipient_hap = 0;
      recipient_hap -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
      recipient_hap = ~recipient_hap; // So H below will now be 1-H
      __m256i _recipient_hap = _mm256_set1_epi32(recipient_hap);
#endif

      f = 0.0; // For accumulating scalar part in ragged end ...
      __m256d _f = _mm256_set1_pd(0.0); // ... and for vector part.

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ alphaRow = &(alpha[N*recipient_alpha]);

#if KALIS_PI == PI_MATRIX
      const double *__restrict__ PiRow = &(Pi[N*recipient]);
      __m256d _rho = _mm256_set1_pd(rho[l-1]);
#elif KALIS_PI == PI_SCALAR
      const double Pirho  = Pi * rho[l-1];
      __m256d _Pirho  = _mm256_set1_pd(Pirho);
#endif

#if KALIS_MU == MU_VECTOR && !defined(KALIS_1STEP)
      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      // ... now update to (1-H) * (1-2*mu) + mu
      const double muTmp1 = 1.0 - 2.0 * mu[l], muTmp2 = mu[l];
      const __m256d _muTmp1 = _mm256_broadcast_sd(&muTmp1), _muTmp2 = _mm256_broadcast_sd(&muTmp2);
#endif

      for(int_fast32_t donoroff=0; donoroff<N/(32*8); ++donoroff) {
#if !defined(KALIS_1STEP)
        // Load next 256 donors and XOR with recipients
        __m256i _HA = _mm256_xor_si256(_recipient_hap, _mm256_load_si256((__m256i*) &(hap_locus[l][donoroff*8])));
        uint32_t *HA = (uint32_t*) &_HA;

        const uint32_t mask = 16843009;
#endif
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

#if KALIS_PI == PI_MATRIX
          __m256d _pi1    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4);
          __m256d _pi2    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 4);
          __m256d _pi3    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 8);
          __m256d _pi4    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 12);

          _pi1            = _mm256_mul_pd(_pi1, _rho);
          _pi2            = _mm256_mul_pd(_pi2, _rho);
          _pi3            = _mm256_mul_pd(_pi3, _rho);
          _pi4            = _mm256_mul_pd(_pi4, _rho);

          _alpha1         = _mm256_fmadd_pd(_alpha1, _omRhoDivF, _pi1); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha2         = _mm256_fmadd_pd(_alpha2, _omRhoDivF, _pi2); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha3         = _mm256_fmadd_pd(_alpha3, _omRhoDivF, _pi3); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha4         = _mm256_fmadd_pd(_alpha4, _omRhoDivF, _pi4); // (Pi*rho + {(1-rho)/f} * alpha)
#elif KALIS_PI == PI_SCALAR
          _alpha1         = _mm256_fmadd_pd(_alpha1, _omRhoDivF, _Pirho); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha2         = _mm256_fmadd_pd(_alpha2, _omRhoDivF, _Pirho); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha3         = _mm256_fmadd_pd(_alpha3, _omRhoDivF, _Pirho); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha4         = _mm256_fmadd_pd(_alpha4, _omRhoDivF, _Pirho); // (Pi*rho + {(1-rho)/f} * alpha)
#endif

#if !defined(KALIS_1STEP)
          __m256d _theta1 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4)/8]) >> (((donor*4)%8)*4), mask))));
          __m256d _theta2 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+1)/8]) >> (((donor*4+1)%8)*4), mask))));
          __m256d _theta3 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+2)/8]) >> (((donor*4+2)%8)*4), mask))));
          __m256d _theta4 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+3)/8]) >> (((donor*4+3)%8)*4), mask))));

          _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#else
          __m256d _theta1 = _mm256_set1_pd(1.0);
          __m256d _theta2 = _mm256_set1_pd(1.0);
          __m256d _theta3 = _mm256_set1_pd(1.0);
          __m256d _theta4 = _mm256_set1_pd(1.0);
#endif

#if KALIS_PI == PI_MATRIX
          _alpha1         = _mm256_mul_pd(_theta1, _alpha1);
          _alpha2         = _mm256_mul_pd(_theta2, _alpha2);
          _alpha3         = _mm256_mul_pd(_theta3, _alpha3);
          _alpha4         = _mm256_mul_pd(_theta4, _alpha4);
#elif KALIS_PI == PI_SCALAR
          _alpha1         = _mm256_and_pd(_mm256_mul_pd(_theta1, _alpha1), _mm256_cmp_pd(_donornum1, _recipient, _CMP_NEQ_UQ));
          _alpha2         = _mm256_and_pd(_mm256_mul_pd(_theta2, _alpha2), _mm256_cmp_pd(_donornum2, _recipient, _CMP_NEQ_UQ));
          _alpha3         = _mm256_and_pd(_mm256_mul_pd(_theta3, _alpha3), _mm256_cmp_pd(_donornum3, _recipient, _CMP_NEQ_UQ));
          _alpha4         = _mm256_and_pd(_mm256_mul_pd(_theta4, _alpha4), _mm256_cmp_pd(_donornum4, _recipient, _CMP_NEQ_UQ));
#endif
          _f              = _mm256_add_pd(_f, _alpha1);
          _f              = _mm256_add_pd(_f, _alpha2);
          _f              = _mm256_add_pd(_f, _alpha3);
          _f              = _mm256_add_pd(_f, _alpha4);

          _mm256_storeu_pd(alphaNow1, _alpha1);
          _mm256_storeu_pd(alphaNow2, _alpha2);
          _mm256_storeu_pd(alphaNow3, _alpha3);
          _mm256_storeu_pd(alphaNow4, _alpha4);

#if KALIS_PI == PI_SCALAR
          _donornum1 = _mm256_add_pd(_donornum1, _donorinc);
          _donornum2 = _mm256_add_pd(_donornum2, _donorinc);
          _donornum3 = _mm256_add_pd(_donornum3, _donorinc);
          _donornum4 = _mm256_add_pd(_donornum4, _donorinc);
#endif
        }
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=0; donor<N%(32*8); ++donor) {
#if !defined(KALIS_1STEP)
        int32_t donor_hap = (hap_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
#endif
        const int32_t donornum = (N/(32*8))*32*8+donor;

#if KALIS_PI == PI_SCALAR && !defined(KALIS_1STEP)
        f += alphaRow[(N/(32*8))*32*8+donor] = (donornum != recipient) * (H * muTmp1 + muTmp2) * (Pirho + omRhoDivF * alphaRow[donornum]);
#elif KALIS_PI == PI_MATRIX && !defined(KALIS_1STEP)
        f += alphaRow[(N/(32*8))*32*8+donor] =                           (H * muTmp1 + muTmp2) * (PiRow[donornum] * rho[l-1] + omRhoDivF * alphaRow[donornum]);
#elif KALIS_PI == PI_SCALAR && defined(KALIS_1STEP)
        f += alphaRow[(N/(32*8))*32*8+donor] = (donornum != recipient) * 1.0 * (Pirho + omRhoDivF * alphaRow[donornum]);
#elif KALIS_PI == PI_MATRIX && defined(KALIS_1STEP)
        f += alphaRow[(N/(32*8))*32*8+donor] =                           1.0 * (PiRow[donornum] * rho[l-1] + omRhoDivF * alphaRow[donornum]);
#endif
      }

      // Accumulate row sum into f
      _f = _mm256_hadd_pd(_f, _f);
      f += ((double*)&_f)[0] + ((double*)&_f)[2];

      fold[recipient_alpha] = f;
    }
  }

#endif

}



void CPP_FN(EXACTFORWARDNOEXP)(NumericMatrix alpha,
            NumericVector alpha_f,
            NumericVector alpha_f2,
            const int alpha_from_rec,
            const int alpha_t,
            const int t,
            const int from_rec,
            const int to_rec,
            const int L,
            const int N,
            PI_TYPE_CPP Pi,
            MU_TYPE_CPP mu,
            NumericVector rho) {
  CPP_RAW_FN(EXACTFORWARDNOEXP)(&(alpha[0]),
             &(alpha_f[0]),
             &(alpha_f2[0]),
             alpha_from_rec,
             alpha_t,
             t,
             from_rec,
             to_rec,
             L,
             N,
             PI_ARG_CPP,
             MU_ARG_CPP,
             &(rho[0]));
}



void PAR_CPP_FN(EXACTFORWARDNOEXP)(NumericMatrix alpha,
                NumericVector alpha_f,
                NumericVector alpha_f2,
                const int alpha_from_rec,
                const int alpha_t,
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
    threads.push_back(std::thread(CPP_RAW_FN(EXACTFORWARDNOEXP),
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
