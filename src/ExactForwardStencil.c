#define _GNU_SOURCE
#include <pthread.h>
#include <stdio.h>
#include <math.h>

#include "Cache.h"



#if defined(KALIS_IMMINTRIN_H)
#include <immintrin.h>
#endif

#if defined(KALIS_ARM_NEON_H)
#include <arm_neon.h>
#endif



#if defined(KALIS_MU) && defined(KALIS_PI)

#include "Stencil2.h"

// Use stencil variables to setup macro arguments for constructing function names
#ifdef KALIS_1STEP
#define OS_FN _1step
#else
#define OS_FN
#endif

#ifdef KALIS_SPEIDEL
#define SP_FN _speidel
#else
#define SP_FN
#endif

#if KALIS_MU == MU_SCALAR
#define MU_FN _scalarMu
#else
#define MU_FN
#endif

#if KALIS_PI == PI_SCALAR
#define PI_FN _scalarPi
#else
#define PI_FN
#endif



struct FWD_CORE_ARGS(OS_FN,SP_FN,MU_FN,PI_FN) {
  double *const restrict alpha;
  double *const restrict alpha_f;
  const size_t alpha_from_rec;
  const size_t alpha_t;
  const size_t t;
  const size_t L;
  const size_t N;
  const PI_TYPE_C Pi;
  const MU_TYPE_C mu;
  const double *const restrict rho;
};
struct FWD_ARGS(OS_FN,SP_FN,MU_FN,PI_FN) {
  struct FWD_CORE_ARGS(OS_FN,SP_FN,MU_FN,PI_FN) *core_args;
  size_t from_rec;
  size_t to_rec;
};

void* FWD_RAW_FN(OS_FN,SP_FN,MU_FN,PI_FN)(void *args) {
  struct FWD_ARGS(OS_FN,SP_FN,MU_FN,PI_FN) *fwd_args;
  fwd_args = (struct FWD_ARGS(OS_FN,SP_FN,MU_FN,PI_FN) *) args;
  double *const restrict alpha = fwd_args->core_args->alpha;
  double *const restrict alpha_f = fwd_args->core_args->alpha_f;
  const size_t alpha_from_rec = fwd_args->core_args->alpha_from_rec;
  const size_t alpha_t = fwd_args->core_args->alpha_t;
  const size_t t = fwd_args->core_args->t;
  const size_t from_rec = fwd_args->from_rec;
  const size_t to_rec = fwd_args->to_rec;
  const size_t L = fwd_args->core_args->L;
  const size_t N = fwd_args->core_args->N;
  const PI_TYPE_C Pi = fwd_args->core_args->Pi;
  const MU_TYPE_C mu = fwd_args->core_args->mu;
  const double *const restrict rho = fwd_args->core_args->rho;

  size_t l = alpha_t;

  double f;
  double *restrict fold;
  fold    = &(alpha_f[0]);

#if KALIS_PI == PI_SCALAR
  double *const restrict PiRow = (double*) malloc(sizeof(double)*N);
  for(size_t i=0; i<N; i++) {
    PiRow[i] = Pi;
  }
#endif

  // Locus zero setup
  if(l>L-1) {
    l = 0;
    for(size_t recipient=from_rec; recipient<to_rec; ++recipient) {
      size_t recipient_alpha = recipient-alpha_from_rec;
      int32_t recipient_hap = (hap_locus[0][recipient/32] >> recipient%32) & 1;
      recipient_hap = 1-recipient_hap; // So H below will now be 1-H

      fold[recipient_alpha]    = 0.0;

      for(size_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (hap_locus[0][donor/32] >> donor%32) & 1;
#ifdef KALIS_SPEIDEL
        int32_t H = (recipient_hap & ~donor_hap) & 1;
#else
        int32_t H = (recipient_hap ^ donor_hap) & 1;
#endif
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

  const size_t reset_l = l;


#if KALIS_MU == MU_SCALAR && !defined(KALIS_1STEP)
  // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
  // ... now update to (1-H) * (1-2*mu) + mu
  const double muTmp1 = 1.0 - 2.0 * mu, muTmp2 = mu;
  const KALIS_DOUBLE _muTmp1 = KALIS_SET_DOUBLE(muTmp1), _muTmp2 = KALIS_SET_DOUBLE(muTmp2);
#endif

  for(size_t recipient=from_rec; recipient<to_rec; ++recipient) {
    size_t recipient_alpha = recipient-alpha_from_rec;
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
      double *restrict alphaRow = &(alpha[N*recipient_alpha]);

#if KALIS_PI == PI_MATRIX
      const double *restrict PiRow = &(Pi[N*recipient]);
#endif
      KALIS_DOUBLE _rho = KALIS_SET_DOUBLE(rho[l-1]);

#if KALIS_MU == MU_VECTOR && !defined(KALIS_1STEP)
      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      // ... now update to (1-H) * (1-2*mu) + mu
      const double muTmp1 = 1.0 - 2.0 * mu[l], muTmp2 = mu[l];
      const KALIS_DOUBLE _muTmp1 = KALIS_SET_DOUBLE(muTmp1), _muTmp2 = KALIS_SET_DOUBLE(muTmp2);
#endif

      for(size_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
#if !defined(KALIS_1STEP)
        // Load next 256 donors and XOR or ANDNOT with recipients
#ifdef KALIS_SPEIDEL
        KALIS_INT32 _HA = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#else
        KALIS_INT32 _HA = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#endif
        int32_t HA[KALIS_INTVEC_SIZE];
        KALIS_STORE_INT_VEC(HA, _HA);
#endif

        // __asm volatile("# LLVM-MCA-BEGIN");
        for(size_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
          // IACA_START
#include KALIS_FORWARD_INNER_UNROLLED(KALIS_UNROLL)
        }
        // __asm volatile("# LLVM-MCA-END");
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(size_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
#if !defined(KALIS_1STEP)
        int32_t donor_hap = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
#ifdef KALIS_SPEIDEL
        int32_t H = (recipient_hap & ~donor_hap) & 1;
#else
        int32_t H = (recipient_hap ^ donor_hap) & 1;
#endif
#endif
        const size_t donornum = (N/(32*KALIS_INTVEC_SIZE))*(32*KALIS_INTVEC_SIZE)+donor;

#if !defined(KALIS_1STEP)
        f += alphaRow[donornum] = (H * muTmp1 + muTmp2) * (PiRow[donornum] * rho[l-1] + omRhoDivF * alphaRow[donornum]);
#elif defined(KALIS_1STEP)
        f += alphaRow[donornum] = 1.0 * (PiRow[donornum] * rho[l-1] + omRhoDivF * alphaRow[donornum]);
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

  return(NULL);
}



void FWD_FN(OS_FN,SP_FN,MU_FN,PI_FN)(double *const restrict alpha,
                               double *const restrict alpha_f,
                               const size_t alpha_from_rec,
                               const size_t alpha_t,
                               const size_t t,
                               const size_t from_rec,
                               const size_t to_rec,
                               const size_t L,
                               const size_t N,
                               const PI_TYPE_C Pi,
                               const MU_TYPE_C mu,
                               const double *const restrict rho,
                               const int *const restrict nthreads,
                               const int nthreads_len) {
  int numthreads, affinity;
  if(nthreads_len > 1) {
    numthreads = nthreads_len;
    affinity = 1;
  } else {
    numthreads = *nthreads;
    affinity = 0;
  }

  struct FWD_CORE_ARGS(OS_FN,SP_FN,MU_FN,PI_FN) fwd_core_args = {
    .alpha = alpha,
    .alpha_f = alpha_f,
    .alpha_from_rec = alpha_from_rec,
    .alpha_t = alpha_t,
    .t = t,
    .L = L,
    .N = N,
    .Pi = Pi,
    .mu = mu,
    .rho = rho
  };
  struct FWD_ARGS(OS_FN,SP_FN,MU_FN,PI_FN) fwd_args[numthreads];
  for(int i=0; i<numthreads; i++) {
    fwd_args[i].core_args = &fwd_core_args;
    fwd_args[i].from_rec = from_rec;
    fwd_args[i].to_rec = to_rec;
  }

  if(numthreads == 1) {
    FWD_RAW_FN(OS_FN,SP_FN,MU_FN,PI_FN)((void*) &fwd_args[0]);
  } else {
    pthread_t threads[numthreads];
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // round(hap(0, nthreads) * double(to_rec-from_rec) / double(nthreads)) + from_rec;
    double spacing = ((double) to_rec-from_rec) / ((double) numthreads);

    for(int i=0; i<numthreads; i++) {
      fwd_args[i].from_rec = (size_t) round(from_rec + i*spacing);
      fwd_args[i].to_rec = (size_t) round(from_rec + (i+1)*spacing);
      pthread_create(&threads[i], &attr, FWD_RAW_FN(OS_FN,SP_FN,MU_FN,PI_FN), (void*) &fwd_args[i]);

#if defined(KALIS_AFFINITY)
      if(affinity) {
        cpu_set_t cpus;
        CPU_ZERO(&cpus);
        CPU_SET(nthreads[i], &cpus);
        int afer = pthread_setaffinity_np(threads[i], sizeof(cpu_set_t), &cpus);
        if(afer != 0)
          printf("Error setting thread affinity!\n");
      }
#endif
    }

    for(int i=0; i<numthreads; i++) {
      pthread_join(threads[i], NULL);
    }
    pthread_attr_destroy(&attr);
  }
}

#endif
