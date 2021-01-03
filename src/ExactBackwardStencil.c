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



struct BCK_CORE_ARGS(SP_FN,MU_FN,PI_FN) {
  double *const restrict beta;
  double *const restrict beta_g;
  const int *const restrict cur_beta_theta;
  const int *const restrict end_beta_theta;
  const size_t beta_from_rec;
  const size_t beta_t;
  const size_t t;
  const size_t L;
  const size_t N;
  PI_TYPE_C Pi;
  const MU_TYPE_C mu;
  const double *const restrict rho;
};
struct BCK_ARGS(SP_FN,MU_FN,PI_FN) {
  struct BCK_CORE_ARGS(SP_FN,MU_FN,PI_FN) *core_args;
  size_t from_rec;
  size_t to_rec;
};

void* BCK_RAW_FN(SP_FN,MU_FN,PI_FN)(void *args) {
  struct BCK_ARGS(SP_FN,MU_FN,PI_FN) *bck_args;
  bck_args = (struct BCK_ARGS(SP_FN,MU_FN,PI_FN) *) args;
  double *const restrict beta = bck_args->core_args->beta;
  double *const restrict beta_g = bck_args->core_args->beta_g;
  const int *const restrict cur_beta_theta = bck_args->core_args->cur_beta_theta;
  const int *const restrict end_beta_theta = bck_args->core_args->end_beta_theta;
  const size_t beta_from_rec = bck_args->core_args->beta_from_rec;
  const size_t beta_t = bck_args->core_args->beta_t;
  const size_t t = bck_args->core_args->t;
  const size_t from_rec = bck_args->from_rec;
  const size_t to_rec = bck_args->to_rec;
  const size_t L = bck_args->core_args->L;
  const size_t N = bck_args->core_args->N;
  PI_TYPE_C Pi = bck_args->core_args->Pi;
  const MU_TYPE_C mu = bck_args->core_args->mu;
  const double *const restrict rho = bck_args->core_args->rho;

  size_t l = beta_t;

  double g;
  double *restrict gold;
  gold    = &(beta_g[0]);

#if KALIS_PI == PI_SCALAR
  double *const restrict PiRow = (double*) malloc(sizeof(double)*N);
  for(size_t i=0; i<N; i++) {
    PiRow[i] = Pi;
  }
#endif

  // Locus L setup
  if(l>L-1) {
    l = L-1;
    for(size_t recipient=from_rec; recipient<to_rec; ++recipient) {
      size_t recipient_beta = recipient-beta_from_rec;
      int32_t recipient_hap = (hap_locus[l][recipient/32] >> recipient%32) & 1;
      recipient_hap = 1-recipient_hap; // So H below will now be 1-H

      gold[recipient_beta] = 0.0;

      for(size_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (hap_locus[l][donor/32] >> donor%32) & 1;
#ifdef KALIS_SPEIDEL
        int32_t H = (recipient_hap & ~donor_hap) & 1;
#else
        int32_t H = (recipient_hap ^ donor_hap) & 1;
#endif
#if KALIS_MU == MU_SCALAR
        double theta = (H * (1.0 - 2.0*mu) + mu);
#elif KALIS_MU == MU_VECTOR
        double theta = (H * (1.0 - 2.0*mu[l]) + mu[l]);
#endif

        beta[donor + N*recipient_beta] = 1.0;
        beta[recipient + N*recipient_beta] = 0.0;

#if KALIS_PI == PI_SCALAR
        gold[recipient_beta] += Pi * theta;
#elif KALIS_PI == PI_MATRIX
        gold[recipient_beta] += Pi[donor + N*recipient] * theta;
#endif
      }

#if KALIS_PI == PI_SCALAR
      // Adjust due to scalar Pi for recipient/receipient locus
      // We can do this at the first locus but not later due to
      // numerical stability problems if (1-mu)*Pi is large
      // compared to all other contributions
#if KALIS_MU == MU_SCALAR
      gold[recipient_beta] -= (1.0-mu)*Pi;
#elif KALIS_MU == MU_VECTOR
      gold[recipient_beta] -= (1.0-mu[0])*Pi;
#endif
#endif
    }
  }

  const size_t reset_l = l;

//#if KALIS_PI == PI_SCALAR
//  const KALIS_DOUBLE _Pi = KALIS_SET_DOUBLE(Pi);
//#endif

#if KALIS_MU == MU_SCALAR
  // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
  const double muTmp1 = 1.0 - 2.0 * mu, muTmp2 = mu;
  const KALIS_DOUBLE _muTmp1 = KALIS_SET_DOUBLE(muTmp1), _muTmp2 = KALIS_SET_DOUBLE(muTmp2);
#endif

  if(l == t && !*cur_beta_theta && *end_beta_theta) {
    // DEBUG: Rcout << "B -> BT (no move, l=" << l << ", all rec)" << std::endl;
    for(size_t recipient=from_rec; recipient<to_rec; ++recipient) {
      size_t recipient_beta = recipient-beta_from_rec;
#if KALIS_PI == PI_SCALAR
      PiRow[recipient] = 0.0;
#endif

      g = 0.0; // For accumulating scalar part in ragged end ...
      KALIS_DOUBLE _g = KALIS_SET_DOUBLE(0.0); // ... and for vector part.

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *const restrict betaRow = &(beta[N*recipient_beta]);
#if KALIS_PI == PI_MATRIX
      double *const restrict PiRow = &(Pi[N*recipient]);
      // The diagonal should be zero ... we can fix that after the loop, but in order to avoid accumulating
      // incorrectly in g, we want to force the diagonal of Pi to zero
      PiRow[recipient] = 0.0;
#endif

      // Load this recipient's bit into all 256-bits of an AVX register
      int32_t recipient_hap = 0;
      recipient_hap -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
      recipient_hap  = ~recipient_hap; // So H below will now be 1-H
      KALIS_INT32 _recipient_hap      = KALIS_SET_INT32(recipient_hap);

#if KALIS_MU == MU_VECTOR
      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      const double muTmp1a = 1.0 - 2.0 * mu[l], muTmp2a = mu[l];
      const KALIS_DOUBLE _muTmp1a = KALIS_SET_DOUBLE(muTmp1a), _muTmp2a = KALIS_SET_DOUBLE(muTmp2a);
#endif

      for(size_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
        // Load next 256 donors and XOR/ANDNOT with recipients
#ifdef KALIS_SPEIDEL
        KALIS_INT32 _HA = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#else
        KALIS_INT32 _HA = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#endif
        int32_t HA[KALIS_INTVEC_SIZE];
        KALIS_STORE_INT_VEC(HA, _HA);

        for(size_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
          // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,E)
        }
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(size_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
        int32_t donor_hapA = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
#ifdef KALIS_SPEIDEL
        int32_t HA = (recipient_hap & ~donor_hapA) & 1;
#else
        int32_t HA = (recipient_hap ^ donor_hapA) & 1;
#endif

        const size_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

#if KALIS_MU == MU_SCALAR
        betaRow[donornum] = (HA * muTmp1 + muTmp2) * betaRow[donornum];
#elif KALIS_MU == MU_VECTOR
        betaRow[donornum] = (HA * muTmp1a + muTmp2a) * betaRow[donornum];
#endif

        g += PiRow[donornum] * betaRow[donornum];
      }

      betaRow[recipient] = 0.0;

      // Accumulate row sum into g
      KALIS_HSUM_DOUBLE(_g);
      g += ((double*)&_g)[0];

      gold[recipient_beta] = g;

#if KALIS_PI == PI_SCALAR
      PiRow[recipient] = Pi;
#endif
    }
  } else if(!(l-1>t) && !*cur_beta_theta && !*end_beta_theta) {
    // If we do one step from entry point, then we have to use the slower code
    // which constructs two thetas
    // ie we can't move into beta*theta space
    //   beta = theta[l]*(beta*theta[l+1] / g...)
    //   g += beta*pi
    // DEBUG: Rcout << "B -> B (l=" << l << ", all r, one step)" << std::endl;
    for(size_t recipient=from_rec; recipient<to_rec; ++recipient) {
      size_t recipient_beta = recipient-beta_from_rec;
#if KALIS_PI == PI_SCALAR
      PiRow[recipient] = 0.0;
#endif
      l = reset_l;

      while(l>t) {
        --l;
        // (1.0-rho[l]) * gratio .... for all recipients we consider to get us going
        double omRhoDivG = (1.0 - rho[l]) / gold[recipient_beta];
        KALIS_DOUBLE _omRhoDivG = KALIS_SET_DOUBLE(omRhoDivG);

        // Load this recipient's bit into all 256-bits of an AVX register
        int32_t recipient_hap = 0, recipient_hap_prev = 0;
        recipient_hap_prev -= (hap_locus[l+1][recipient/32] >> recipient%32) & 1;
        recipient_hap_prev  = ~recipient_hap_prev; // So H below will now be 1-H
        recipient_hap      -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
        recipient_hap       = ~recipient_hap; // So H below will now be 1-H
        KALIS_INT32 _recipient_hap_prev = KALIS_SET_INT32(recipient_hap_prev);
        KALIS_INT32 _recipient_hap      = KALIS_SET_INT32(recipient_hap);

        g = 0.0; // For accumulating scalar part in ragged end ...
        KALIS_DOUBLE _g = KALIS_SET_DOUBLE(0.0); // ... and for vector part.

        // TODO: for larger problems break this down into L1 cachable chunks of
        //       donors at a time
        double *const restrict betaRow = &(beta[N*recipient_beta]);
#if KALIS_PI == PI_MATRIX
        double *const restrict PiRow = &(Pi[N*recipient]);
        // The diagonal should be zero ... we can fix that after the loop, but in order to avoid accumulating
        // incorrectly in g, we want to force the diagonal of Pi to zero
        PiRow[recipient] = 0.0;
#endif

#if KALIS_MU == MU_VECTOR
        // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
        const double muTmp1a = 1.0 - 2.0 * mu[l+1], muTmp2a = mu[l+1];
        const KALIS_DOUBLE _muTmp1a = KALIS_SET_DOUBLE(muTmp1a), _muTmp2a = KALIS_SET_DOUBLE(muTmp2a);
        const double muTmp1b = 1.0 - 2.0 * mu[l], muTmp2b = mu[l];
        const KALIS_DOUBLE _muTmp1b = KALIS_SET_DOUBLE(muTmp1b), _muTmp2b = KALIS_SET_DOUBLE(muTmp2b);
#endif

        // Setup rho for AVX ops
        KALIS_DOUBLE _rho = KALIS_SET_DOUBLE(rho[l]);
        for(size_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
          // Load next 256 donors and XOR/ANDNOT with recipients
#ifdef KALIS_SPEIDEL
          KALIS_INT32 _HA = KALIS_ANDNOT_INT(_recipient_hap_prev, KALIS_LOAD_INT_VEC(hap_locus[l+1][donoroff*KALIS_INTVEC_SIZE]));
          KALIS_INT32 _HB = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#else
          KALIS_INT32 _HA = KALIS_XOR_INT(_recipient_hap_prev, KALIS_LOAD_INT_VEC(hap_locus[l+1][donoroff*KALIS_INTVEC_SIZE]));
          KALIS_INT32 _HB = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#endif
          int32_t HA[KALIS_INTVEC_SIZE], HB[KALIS_INTVEC_SIZE];
          KALIS_STORE_INT_VEC(HA, _HA);
          KALIS_STORE_INT_VEC(HB, _HB);

          for(size_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
            // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,A)
          }
          // IACA_END
        }
        // Tidy up any ragged end past a multiple of 256 ...
        for(size_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
          int32_t donor_hap = (hap_locus[l+1][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
#ifdef KALIS_SPEIDEL
          int32_t H = (recipient_hap_prev & ~donor_hap) & 1;
#else
          int32_t H = (recipient_hap_prev ^ donor_hap) & 1;
#endif

          const size_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

#if KALIS_MU == MU_SCALAR
          betaRow[donornum] = rho[l] + (H * muTmp1 + muTmp2) * betaRow[donornum] * omRhoDivG;
#elif KALIS_MU == MU_VECTOR
          betaRow[donornum] = rho[l] + (H * muTmp1a + muTmp2a) * betaRow[donornum] * omRhoDivG;
#endif

          donor_hap = (hap_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
#ifdef KALIS_SPEIDEL
          H = (recipient_hap & ~donor_hap) & 1;
#else
          H = (recipient_hap ^ donor_hap) & 1;
#endif

#if KALIS_MU == MU_VECTOR
          g += PiRow[donornum] * (H * muTmp1b + muTmp2b) * betaRow[donornum];
#elif KALIS_MU == MU_SCALAR
          g += PiRow[donornum] * (H * muTmp1 + muTmp2) * betaRow[donornum];
#endif
        }
        betaRow[recipient] = 0.0;

        // Accumulate row sum into g
        KALIS_HSUM_DOUBLE(_g);
        g += ((double*)&_g)[0];

        gold[recipient_beta] = g;
      }

#if KALIS_PI == PI_SCALAR
      PiRow[recipient] = Pi;
#endif
    }
  } else {
    // If we do more than one step from entry point, then for the intermediate
    // steps we track theta*beta to avoid double theta computation
    // ie on first step, move into theta*beta space at the same time as the
    //   locus step; then march back staying in theta*beta space; then on
    //   final step move back to beta space at the same time as the locus step
    // See if clauses below for detail
    for(size_t recipient=from_rec; recipient<to_rec; ++recipient) {
      size_t recipient_beta = recipient-beta_from_rec;
#if KALIS_PI == PI_SCALAR
      PiRow[recipient] = 0.0;
#endif
      l = reset_l;

      while(l>t) {
        --l;
        // (1.0-rho[l]) * gratio .... for all recipients we consider to get us going
        double omRhoDivG = (1.0 - rho[l]) / gold[recipient_beta];
        KALIS_DOUBLE _omRhoDivG = KALIS_SET_DOUBLE(omRhoDivG);

        g = 0.0; // For accumulating scalar part in ragged end ...
        KALIS_DOUBLE _g = KALIS_SET_DOUBLE(0.0); // ... and for vector part.

        // TODO: for larger problems break this down into L1 cachable chunks of
        //       donors at a time
        double *const restrict betaRow = &(beta[N*recipient_beta]);
#if KALIS_PI == PI_MATRIX
        double *const restrict PiRow = &(Pi[N*recipient]);
        // The diagonal should be zero ... we can fix that after the loop, but in order to avoid accumulating
        // incorrectly in g, we want to force the diagonal of Pi to zero
        PiRow[recipient] = 0.0;
#endif

        // Setup rho for AVX ops
        KALIS_DOUBLE _rho = KALIS_SET_DOUBLE(rho[l]);

        if((l<reset_l-1 || *cur_beta_theta) && (l>t || *end_beta_theta)) {
          // Mid beta*theta step
          // Perform locus updates staying in beta*theta space
          //   beta = theta[l-i]*(beta / g...)
          //   g += beta*pi
          // DEBUG: Rcout << "BT -> BT (l=" << l <<", r=" << recipient << ")" << std::endl;

          // Load this recipient's bit into all 256-bits of an AVX register
          int32_t recipient_hap = 0;
          recipient_hap -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
          recipient_hap  = ~recipient_hap; // So H below will now be 1-H
          KALIS_INT32 _recipient_hap = KALIS_SET_INT32(recipient_hap);

#if KALIS_MU == MU_VECTOR
          // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
          const double muTmp1b = 1.0 - 2.0 * mu[l], muTmp2b = mu[l];
          const KALIS_DOUBLE _muTmp1b = KALIS_SET_DOUBLE(muTmp1b), _muTmp2b = KALIS_SET_DOUBLE(muTmp2b);
#endif

          for(size_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
            // Load next 256 donors and XOR/ANDNOT with recipients
#ifdef KALIS_SPEIDEL
            KALIS_INT32 _HB = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#else
            KALIS_INT32 _HB = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#endif
            int32_t HB[KALIS_INTVEC_SIZE];
            KALIS_STORE_INT_VEC(HB, _HB);

            for(size_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
              // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,B)
            }
            // IACA_END
          }
          // Tidy up any ragged end past a multiple of 256 ...
          for(size_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
            int32_t donor_hap = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
#ifdef KALIS_SPEIDEL
            int32_t H = (recipient_hap & ~donor_hap) & 1;
#else
            int32_t H = (recipient_hap ^ donor_hap) & 1;
#endif

            const size_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

#if KALIS_MU == MU_SCALAR
            betaRow[donornum] = (rho[l] + betaRow[donornum] * omRhoDivG) * (H * muTmp1 + muTmp2);
#elif KALIS_MU == MU_VECTOR
            betaRow[donornum] = (rho[l] + betaRow[donornum] * omRhoDivG) * (H * muTmp1b + muTmp2b);
#endif

            g += PiRow[donornum] * betaRow[donornum];
          }

          betaRow[recipient] = 0.0;

        } else if(l==reset_l-1 && !*cur_beta_theta) {
          // Pre beta*theta step
          // Move from beta space into beta*theta space simultaneously with a locus update
          //   beta = theta[l-1]*(beta / g...)
          //   g += beta*pi
          // DEBUG: Rcout << "B -> BT (l=" << l <<", r=" << recipient << ")" << std::endl;

          // Load this recipient's bit into all 256-bits of an AVX register
          int32_t recipient_hap = 0, recipient_hap_prev = 0;
          recipient_hap_prev -= (hap_locus[l+1][recipient/32] >> recipient%32) & 1;
          recipient_hap_prev  = ~recipient_hap_prev; // So H below will now be 1-H
          recipient_hap      -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
          recipient_hap       = ~recipient_hap; // So H below will now be 1-H
          KALIS_INT32 _recipient_hap_prev = KALIS_SET_INT32(recipient_hap_prev);
          KALIS_INT32 _recipient_hap      = KALIS_SET_INT32(recipient_hap);

#if KALIS_MU == MU_VECTOR
          // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
          const double muTmp1a = 1.0 - 2.0 * mu[l+1], muTmp2a = mu[l+1];
          const KALIS_DOUBLE _muTmp1a = KALIS_SET_DOUBLE(muTmp1a), _muTmp2a = KALIS_SET_DOUBLE(muTmp2a);
          const double muTmp1b = 1.0 - 2.0 * mu[l], muTmp2b = mu[l];
          const KALIS_DOUBLE _muTmp1b = KALIS_SET_DOUBLE(muTmp1b), _muTmp2b = KALIS_SET_DOUBLE(muTmp2b);
#endif

          for(size_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
            // Load next 256 donors and XOR/ANDNOT with recipients
#ifdef KALIS_SPEIDEL
            KALIS_INT32 _HA = KALIS_ANDNOT_INT(_recipient_hap_prev, KALIS_LOAD_INT_VEC(hap_locus[l+1][donoroff*KALIS_INTVEC_SIZE]));
            KALIS_INT32 _HB = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#else
            KALIS_INT32 _HA = KALIS_XOR_INT(_recipient_hap_prev, KALIS_LOAD_INT_VEC(hap_locus[l+1][donoroff*KALIS_INTVEC_SIZE]));
            KALIS_INT32 _HB = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#endif
            int32_t HA[KALIS_INTVEC_SIZE], HB[KALIS_INTVEC_SIZE];
            KALIS_STORE_INT_VEC(HA, _HA);
            KALIS_STORE_INT_VEC(HB, _HB);

            for(size_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
              // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,C)
            }
            // IACA_END
          }
          // Tidy up any ragged end past a multiple of 256 ...
          for(size_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
            int32_t donor_hapA = (hap_locus[l+1][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
            int32_t donor_hap = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
#ifdef KALIS_SPEIDEL
            int32_t HA = (recipient_hap_prev & ~donor_hapA) & 1;
            int32_t H  = (recipient_hap & ~donor_hap) & 1;
#else
            int32_t HA = (recipient_hap_prev ^ donor_hapA) & 1;
            int32_t H  = (recipient_hap ^ donor_hap) & 1;
#endif

            const size_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

#if KALIS_MU == MU_SCALAR
            betaRow[donornum] = (rho[l] + (HA * muTmp1 + muTmp2) * betaRow[donornum] * omRhoDivG) * (H * muTmp1 + muTmp2);
#elif KALIS_MU == MU_VECTOR
            betaRow[donornum] = (rho[l] + (HA * muTmp1a + muTmp2a) * betaRow[donornum] * omRhoDivG) * (H * muTmp1b + muTmp2b);
#endif

            g += PiRow[donornum] * betaRow[donornum];
          }

          betaRow[recipient] = 0.0;

        } else {
          // Exiting beta*theta step
          // Move from beta*theta space into beta space simultaneously with a locus update
          //   beta = beta / g...
          //   g += beta*theta[l-targ]*pi
          // DEBUG: Rcout << "BT -> B (l=" << l <<", r=" << recipient << ")" << std::endl;

          // Load this recipient's bit into all 256-bits of an AVX register
          int32_t recipient_hap = 0;
          recipient_hap -= (hap_locus[l][recipient/32] >> recipient%32) & 1;
          recipient_hap  = ~recipient_hap; // So H below will now be 1-H
          KALIS_INT32 _recipient_hap = KALIS_SET_INT32(recipient_hap);

#if KALIS_MU == MU_VECTOR
          // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
          const double muTmp1b = 1.0 - 2.0 * mu[l], muTmp2b = mu[l];
          const KALIS_DOUBLE _muTmp1b = KALIS_SET_DOUBLE(muTmp1b), _muTmp2b = KALIS_SET_DOUBLE(muTmp2b);
#endif

          for(size_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
            // Load next 256 donors and XOR/ANDNOT with recipients
#ifdef KALIS_SPEIDEL
            KALIS_INT32 _HB = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#else
            KALIS_INT32 _HB = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
#endif
            int32_t HB[KALIS_INTVEC_SIZE];
            KALIS_STORE_INT_VEC(HB, _HB);

            for(size_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
              // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,D)
            }
            // IACA_END
          }
          // Tidy up any ragged end past a multiple of 256 ...
          for(size_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
            int32_t donor_hap = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
#ifdef KALIS_SPEIDEL
            int32_t H = (recipient_hap & ~donor_hap) & 1;
#else
            int32_t H = (recipient_hap ^ donor_hap) & 1;
#endif

            const size_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

#if KALIS_MU == MU_SCALAR
            betaRow[donornum] = (rho[l] + betaRow[donornum] * omRhoDivG);
#elif KALIS_MU == MU_VECTOR
            betaRow[donornum] = (rho[l] + betaRow[donornum] * omRhoDivG);
#endif

#if KALIS_MU == MU_VECTOR
            g += PiRow[donornum] * betaRow[donornum] * (H * muTmp1b + muTmp2b);
#elif KALIS_MU == MU_SCALAR
            g += PiRow[donornum] * betaRow[donornum] * (H * muTmp1 + muTmp2);
#endif
          }

          betaRow[recipient] = 0.0;

        }

        // Accumulate row sum into g
        KALIS_HSUM_DOUBLE(_g);
        g += ((double*)&_g)[0];

        gold[recipient_beta] = g;
      }

#if KALIS_PI == PI_SCALAR
      PiRow[recipient] = Pi;
#endif
    }
  }

#if KALIS_PI == PI_SCALAR
  free(PiRow);
#endif

  return(NULL);
}



void BCK_FN(SP_FN,MU_FN,PI_FN)(double *const restrict beta,
                               double *const restrict beta_g,
                               int *const restrict cur_beta_theta,
                               int *const restrict end_beta_theta,
                               const size_t beta_from_rec,
                               const size_t beta_t,
                               const size_t t,
                               const size_t from_rec,
                               const size_t to_rec,
                               const size_t L,
                               const size_t N,
                               PI_TYPE_C Pi,
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

  struct BCK_CORE_ARGS(SP_FN,MU_FN,PI_FN) bck_core_args = {
    .beta = beta,
    .beta_g = beta_g,
    .cur_beta_theta = cur_beta_theta,
    .end_beta_theta = end_beta_theta,
    .beta_from_rec = beta_from_rec,
    .beta_t = beta_t,
    .t = t,
    .L = L,
    .N = N,
    .Pi = Pi,
    .mu = mu,
    .rho = rho
  };

  struct BCK_ARGS(SP_FN,MU_FN,PI_FN) bck_args[numthreads];
  for(int i=0; i<numthreads; i++) {
    bck_args[i].core_args = &bck_core_args;
    bck_args[i].from_rec = from_rec;
    bck_args[i].to_rec = to_rec;
  }

  if(numthreads == 1) {
    BCK_RAW_FN(SP_FN,MU_FN,PI_FN)((void*) &bck_args[0]);
  } else {
    pthread_t threads[numthreads];
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // round(hap(0, nthreads) * double(to_rec-from_rec) / double(nthreads)) + from_rec;
    double spacing = ((double) to_rec-from_rec) / ((double) numthreads);

    for(int i=0; i<numthreads; i++) {
      bck_args[i].from_rec = (size_t) round(from_rec + i*spacing);
      bck_args[i].to_rec = (size_t) round(from_rec + (i+1)*spacing);
      pthread_create(&threads[i], &attr, BCK_RAW_FN(SP_FN,MU_FN,PI_FN), (void*) &bck_args[i]);

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

  *(&(cur_beta_theta[0])) = *(&(end_beta_theta[0]));
}

#endif
