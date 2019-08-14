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
  int_fast32_t l = alpha_t;

  double f;
  double *__restrict__ fold;
  fold    = &(alpha_f[0]);

#if KALIS_PI == PI_SCALAR
  double *__restrict__ const PiRow = (double*) malloc(sizeof(double)*N);
  for(int_fast32_t i=0; i<N; i++) {
    PiRow[i] = Pi;
  }
#endif

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
  const KALIS_DOUBLE _muTmp1 = KALIS_SET_DOUBLE(muTmp1), _muTmp2 = KALIS_SET_DOUBLE(muTmp2);
#endif

  for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
    int_fast32_t recipient_alpha = recipient-alpha_from_rec;
#if KALIS_PI == PI_SCALAR
    PiRow[recipient] = 0.0;
#endif
    l = reset_l;

    while(l<t) {
      ++l;
      // (1.0-rho[l-1]) * fratio .... for all recipients we consider to get us going
      double omRhoDivF = (1.0 - rho[l-1]) / fold[recipient_alpha];
      KALIS_DOUBLE _omRhoDivF = KALIS_SET_DOUBLE(omRhoDivF);

#if !defined(KALIS_1STEP)
      // Load this recipient's bit into all 256-bits of an AVX register
      int32_t recipient_hap = 0;
      recipient_hap -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
      recipient_hap = ~recipient_hap; // So H below will now be 1-H
      KALIS_INT32 _recipient_hap = KALIS_SET_INT32(recipient_hap);
#endif

      f = 0.0; // For accumulating scalar part in ragged end ...
      KALIS_DOUBLE _f = KALIS_SET_DOUBLE(0.0); // ... and for vector part.

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ alphaRow = &(alpha[N*recipient_alpha]);

#if KALIS_PI == PI_MATRIX
      const double *__restrict__ PiRow = &(Pi[N*recipient]);
#endif
      KALIS_DOUBLE _rho = KALIS_SET_DOUBLE(rho[l-1]);

#if KALIS_MU == MU_VECTOR && !defined(KALIS_1STEP)
      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      // ... now update to (1-H) * (1-2*mu) + mu
      const double muTmp1 = 1.0 - 2.0 * mu[l], muTmp2 = mu[l];
      const KALIS_DOUBLE _muTmp1 = KALIS_SET_DOUBLE(muTmp1), _muTmp2 = KALIS_SET_DOUBLE(muTmp2);
#endif

      for(int_fast32_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
#if !defined(KALIS_1STEP)
        // Load next 256 donors and XOR with recipients
        KALIS_INT32 _HA = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
        uint32_t *HA = (uint32_t*) &_HA;
#endif

        for(int_fast32_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/4; ++donor) {
          // IACA_START
          double *alphaNow1    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*4);
          double *alphaNow2    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*4) + (  KALIS_DOUBLEVEC_SIZE);
          double *alphaNow3    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*4) + (2*KALIS_DOUBLEVEC_SIZE);
          double *alphaNow4    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*4) + (3*KALIS_DOUBLEVEC_SIZE);

          KALIS_DOUBLE _alpha1 = KALIS_LOADU_DOUBLE(alphaNow1);
          KALIS_DOUBLE _alpha2 = KALIS_LOADU_DOUBLE(alphaNow2);
          KALIS_DOUBLE _alpha3 = KALIS_LOADU_DOUBLE(alphaNow3);
          KALIS_DOUBLE _alpha4 = KALIS_LOADU_DOUBLE(alphaNow4);

          KALIS_DOUBLE _pi1    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*4)                           );
          KALIS_DOUBLE _pi2    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*4) + (  KALIS_DOUBLEVEC_SIZE));
          KALIS_DOUBLE _pi3    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*4) + (2*KALIS_DOUBLEVEC_SIZE));
          KALIS_DOUBLE _pi4    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*4) + (3*KALIS_DOUBLEVEC_SIZE));

          _pi1                 = KALIS_MUL_DOUBLE(_pi1, _rho);
          _pi2                 = KALIS_MUL_DOUBLE(_pi2, _rho);
          _pi3                 = KALIS_MUL_DOUBLE(_pi3, _rho);
          _pi4                 = KALIS_MUL_DOUBLE(_pi4, _rho);

          _alpha1              = KALIS_FMA_DOUBLE(_alpha1, _omRhoDivF, _pi1); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha2              = KALIS_FMA_DOUBLE(_alpha2, _omRhoDivF, _pi2); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha3              = KALIS_FMA_DOUBLE(_alpha3, _omRhoDivF, _pi3); // (Pi*rho + {(1-rho)/f} * alpha)
          _alpha4              = KALIS_FMA_DOUBLE(_alpha4, _omRhoDivF, _pi4); // (Pi*rho + {(1-rho)/f} * alpha)

#if !defined(KALIS_1STEP)
          KALIS_DOUBLE _theta1 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*4))/32])                          >> ((donor*(KALIS_DOUBLEVEC_SIZE*4))%32));
          KALIS_DOUBLE _theta2 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*4)+(  KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*4)+(  KALIS_DOUBLEVEC_SIZE))%32));
          KALIS_DOUBLE _theta3 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*4)+(2*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*4)+(2*KALIS_DOUBLEVEC_SIZE))%32));
          KALIS_DOUBLE _theta4 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*4)+(3*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*4)+(3*KALIS_DOUBLEVEC_SIZE))%32));

          _theta1              = KALIS_FMA_DOUBLE(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta2              = KALIS_FMA_DOUBLE(_theta2, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta3              = KALIS_FMA_DOUBLE(_theta3, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
          _theta4              = KALIS_FMA_DOUBLE(_theta4, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#else
          KALIS_DOUBLE _theta1 = KALIS_SET_DOUBLE(1.0);
          KALIS_DOUBLE _theta2 = KALIS_SET_DOUBLE(1.0);
          KALIS_DOUBLE _theta3 = KALIS_SET_DOUBLE(1.0);
          KALIS_DOUBLE _theta4 = KALIS_SET_DOUBLE(1.0);
#endif

          _alpha1              = KALIS_MUL_DOUBLE(_theta1, _alpha1);
          _alpha2              = KALIS_MUL_DOUBLE(_theta2, _alpha2);
          _alpha3              = KALIS_MUL_DOUBLE(_theta3, _alpha3);
          _alpha4              = KALIS_MUL_DOUBLE(_theta4, _alpha4);

          _f                   = KALIS_ADD_DOUBLE(_f, _alpha1);
          _f                   = KALIS_ADD_DOUBLE(_f, _alpha2);
          _f                   = KALIS_ADD_DOUBLE(_f, _alpha3);
          _f                   = KALIS_ADD_DOUBLE(_f, _alpha4);

          KALIS_STOREU_DOUBLE(alphaNow1, _alpha1);
          KALIS_STOREU_DOUBLE(alphaNow2, _alpha2);
          KALIS_STOREU_DOUBLE(alphaNow3, _alpha3);
          KALIS_STOREU_DOUBLE(alphaNow4, _alpha4);
        }
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
#if !defined(KALIS_1STEP)
        int32_t donor_hap = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
#endif
        const int32_t donornum = (N/(32*KALIS_INTVEC_SIZE))*(32*KALIS_INTVEC_SIZE)+donor;

#if !defined(KALIS_1STEP)
        f += alphaRow[(N/(32*KALIS_INTVEC_SIZE))*(32*KALIS_INTVEC_SIZE)+donor] = (H * muTmp1 + muTmp2) * (PiRow[donornum] * rho[l-1] + omRhoDivF * alphaRow[donornum]);
#elif defined(KALIS_1STEP)
        f += alphaRow[(N/(32*KALIS_INTVEC_SIZE))*(32*KALIS_INTVEC_SIZE)+donor] = 1.0 * (PiRow[donornum] * rho[l-1] + omRhoDivF * alphaRow[donornum]);
#endif
      }

      // Accumulate row sum into f
      KALIS_HSUM_DOUBLE(_f);
      f += ((double*)&_f)[0];
      //_f = _mm256_hadd_pd(_f, _f);
      //f += ((double*)&_f)[0] + ((double*)&_f)[2];

      fold[recipient_alpha] = f;
    }

#if KALIS_PI == PI_SCALAR
    PiRow[recipient] = Pi;
#endif
  }

#if KALIS_PI == PI_SCALAR
  free(PiRow);
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
