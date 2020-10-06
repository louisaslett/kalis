#include <Rcpp.h>
using namespace Rcpp;
#include <stdlib.h>
#include <thread>
#include <vector>
#include <functional>

#include "Cache.h"


#if defined(KALIS_IMMINTRIN_H)
#include <immintrin.h>
#endif
#if defined(KALIS_ARM_NEON_H)
#include <arm_neon.h>
#endif

// #include <iacaMarks.h>

#ifdef EXACTBACKWARDNOEXP

#include "Stencil2.h"

void CPP_RAW_FN(EXACTBACKWARDNOEXP)(double *const __restrict__ beta,
                double *const __restrict__ beta_g,
                const int *const __restrict__ cur_beta_theta,
                const int *const __restrict__ end_beta_theta,
                const int beta_from_rec,
                const int beta_t,
                const int t,
                const int from_rec,
                const int to_rec,
                const int L,
                const int N,
                PI_TYPE_C Pi,
                const MU_TYPE_C mu,
                const double *const __restrict__ rho,
                const bool use_speidel) {
  int_fast32_t l = beta_t;
  // DEBUG: Rcout << "Backward called" << std::endl;
  double g;
  double *__restrict__ gold;
  gold    = &(beta_g[0]);

#if KALIS_PI == PI_SCALAR
  double *__restrict__ const PiRow = (double*) malloc(sizeof(double)*N);
  for(int_fast32_t i=0; i<N; i++) {
    PiRow[i] = Pi;
  }
#endif

  // Locus L setup
  if(l>L-1) {
    l = L-1;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
      int32_t recipient_hap = (hap_locus[l][recipient/32] >> recipient%32) & 1;
      recipient_hap = 1-recipient_hap; // So H below will now be 1-H

      gold[recipient_beta] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (hap_locus[l][donor/32] >> donor%32) & 1;
        int32_t H;
        if(use_speidel)
          H = (recipient_hap & ~donor_hap) & 1;
        else
          H = (recipient_hap ^ donor_hap) & 1;
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

  const int_fast32_t reset_l = l;

#if KALIS_PI == PI_SCALAR
  const KALIS_DOUBLE _Pi = KALIS_SET_DOUBLE(Pi);
#endif

#if KALIS_MU == MU_SCALAR
  // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
  const double muTmp1 = 1.0 - 2.0 * mu, muTmp2 = mu;
  const KALIS_DOUBLE _muTmp1 = KALIS_SET_DOUBLE(muTmp1), _muTmp2 = KALIS_SET_DOUBLE(muTmp2);
#endif

  if(l == t && !*cur_beta_theta && *end_beta_theta) {
    // DEBUG: Rcout << "B -> BT (no move, l=" << l << ", all rec)" << std::endl;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
#if KALIS_PI == PI_SCALAR
      PiRow[recipient] = 0.0;
#endif

      g = 0.0; // For accumulating scalar part in ragged end ...
      KALIS_DOUBLE _g = KALIS_SET_DOUBLE(0.0); // ... and for vector part.

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ const betaRow = &(beta[N*recipient_beta]);
#if KALIS_PI == PI_MATRIX
      double *__restrict__ const PiRow = &(Pi[N*recipient]);
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

      for(int_fast32_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
        // Load next 256 donors and XOR/ANDNOT with recipients
        KALIS_INT32 _HA;
        if(use_speidel)
          _HA = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
        else
          _HA = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
        uint32_t *HA = (uint32_t*) &_HA;

        for(int_fast32_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
          // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,E)
        }
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
        int32_t donor_hapA = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
        int32_t HA;
        if(use_speidel)
          HA = (recipient_hap & ~donor_hapA) & 1;
        else
          HA = (recipient_hap ^ donor_hapA) & 1;

        const int32_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

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
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
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
        double *__restrict__ const betaRow = &(beta[N*recipient_beta]);
#if KALIS_PI == PI_MATRIX
        double *__restrict__ const PiRow = &(Pi[N*recipient]);
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
        for(int_fast32_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
          // Load next 256 donors and XOR/ANDNOT with recipients
          KALIS_INT32 _HA, _HB;
          if(use_speidel) {
            _HA = KALIS_ANDNOT_INT(_recipient_hap_prev, KALIS_LOAD_INT_VEC(hap_locus[l+1][donoroff*KALIS_INTVEC_SIZE]));
            _HB = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
          } else {
            _HA = KALIS_XOR_INT(_recipient_hap_prev, KALIS_LOAD_INT_VEC(hap_locus[l+1][donoroff*KALIS_INTVEC_SIZE]));
            _HB = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
          }
          uint32_t *HA = (uint32_t*) &_HA;
          uint32_t *HB = (uint32_t*) &_HB;

          for(int_fast32_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
            // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,A)
          }
          // IACA_END
        }
        // Tidy up any ragged end past a multiple of 256 ...
        for(int32_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
          int32_t donor_hap = (hap_locus[l+1][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
          int32_t H;
          if(use_speidel)
            H = (recipient_hap_prev & ~donor_hap) & 1;
          else
            H = (recipient_hap_prev ^ donor_hap) & 1;

          const int32_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

#if KALIS_MU == MU_SCALAR
          betaRow[donornum] = rho[l] + (H * muTmp1 + muTmp2) * betaRow[donornum] * omRhoDivG;
#elif KALIS_MU == MU_VECTOR
          betaRow[donornum] = rho[l] + (H * muTmp1a + muTmp2a) * betaRow[donornum] * omRhoDivG;
#endif

          donor_hap = (hap_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
          if(use_speidel)
            H = (recipient_hap & ~donor_hap) & 1;
          else
            H = (recipient_hap ^ donor_hap) & 1;

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
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
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
        double *__restrict__ const betaRow = &(beta[N*recipient_beta]);
#if KALIS_PI == PI_MATRIX
        double *__restrict__ const PiRow = &(Pi[N*recipient]);
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

          for(int_fast32_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
            // Load next 256 donors and XOR/ANDNOT with recipients
            KALIS_INT32 _HB;
            if(use_speidel)
              _HB = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
            else
              _HB = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
            uint32_t *HB = (uint32_t*) &_HB;

            for(int_fast32_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
              // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,B)
            }
            // IACA_END
          }
          // Tidy up any ragged end past a multiple of 256 ...
          for(int32_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
            int32_t donor_hap = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
            int32_t H;
            if(use_speidel)
              H = (recipient_hap & ~donor_hap) & 1;
            else
              H = (recipient_hap ^ donor_hap) & 1;

            const int32_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

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

          for(int_fast32_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
            // Load next 256 donors and XOR/ANDNOT with recipients
            KALIS_INT32 _HA, _HB;
            if(use_speidel) {
              _HA = KALIS_ANDNOT_INT(_recipient_hap_prev, KALIS_LOAD_INT_VEC(hap_locus[l+1][donoroff*KALIS_INTVEC_SIZE]));
              _HB = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
            } else {
              _HA = KALIS_XOR_INT(_recipient_hap_prev, KALIS_LOAD_INT_VEC(hap_locus[l+1][donoroff*KALIS_INTVEC_SIZE]));
              _HB = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
            }
            uint32_t *HA = (uint32_t*) &_HA;
            uint32_t *HB = (uint32_t*) &_HB;

            for(int_fast32_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
              // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,C)
            }
            // IACA_END
          }
          // Tidy up any ragged end past a multiple of 256 ...
          for(int32_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
            int32_t donor_hapA = (hap_locus[l+1][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
            int32_t donor_hap = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
            int32_t HA, H;
            if(use_speidel) {
              HA = (recipient_hap_prev & ~donor_hapA) & 1;
              H  = (recipient_hap & ~donor_hap) & 1;
            } else {
              HA = (recipient_hap_prev ^ donor_hapA) & 1;
              H  = (recipient_hap ^ donor_hap) & 1;
            }

            const int32_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

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

          for(int_fast32_t donoroff=0; donoroff<N/(32*KALIS_INTVEC_SIZE); ++donoroff) {
            // Load next 256 donors and XOR/ANDNOT with recipients
            KALIS_INT32 _HB;
            if(use_speidel)
              _HB = KALIS_ANDNOT_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
            else
              _HB = KALIS_XOR_INT(_recipient_hap, KALIS_LOAD_INT_VEC(hap_locus[l][donoroff*KALIS_INTVEC_SIZE]));
            uint32_t *HB = (uint32_t*) &_HB;

            for(int_fast32_t donor=0; donor<((32*KALIS_INTVEC_SIZE)/KALIS_DOUBLEVEC_SIZE)/KALIS_UNROLL; ++donor) {
              // IACA_START
#include KALIS_BACKWARD_INNER_UNROLLED(KALIS_UNROLL,D)
            }
            // IACA_END
          }
          // Tidy up any ragged end past a multiple of 256 ...
          for(int32_t donor=0; donor<N%(32*KALIS_INTVEC_SIZE); ++donor) {
            int32_t donor_hap = (hap_locus[l][(N/(32*KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE + donor/32] >> (donor%32)) & 1;
            int32_t H;
            if(use_speidel)
              H = (recipient_hap & ~donor_hap) & 1;
            else
              H = (recipient_hap ^ donor_hap) & 1;

            const int32_t donornum = (N/(32*KALIS_INTVEC_SIZE))*32*KALIS_INTVEC_SIZE+donor;

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
}



void CPP_FN(EXACTBACKWARDNOEXP)(NumericMatrix beta,
            NumericVector beta_g,
            LogicalVector cur_beta_theta,
            LogicalVector end_beta_theta,
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
            const bool use_speidel) {
  CPP_RAW_FN(EXACTBACKWARDNOEXP)(&(beta[0]),
             &(beta_g[0]),
             &(cur_beta_theta[0]),
             &(end_beta_theta[0]),
             beta_from_rec,
             beta_t,
             t,
             from_rec,
             to_rec,
             L,
             N,
             PI_ARG_CPP,
             MU_ARG_CPP,
             &(rho[0]),
             use_speidel);
  *(&(cur_beta_theta[0])) = *(&(end_beta_theta[0]));
}



void PAR_CPP_FN(EXACTBACKWARDNOEXP)(NumericMatrix beta,
                NumericVector beta_g,
                LogicalVector cur_beta_theta,
                LogicalVector end_beta_theta,
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
                const bool use_speidel,
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
                                  &(cur_beta_theta[0]),
                                  &(end_beta_theta[0]),
                                  beta_from_rec,
                                  beta_t,
                                  t,
                                  round(from_rec + i*spacing),
                                  round(from_rec + (i+1)*spacing),
                                  L,
                                  N,
                                  PI_ARG_CPP,
                                  MU_ARG_CPP,
                                  &(rho[0]),
                                  use_speidel));
    // Rcout << "From: " << round(from_rec + i*spacing) << ", To: " << round(from_rec + (i+1)*spacing) << "\n";
  }

  for(auto& th : threads) {
    th.join();
  }

  *(&(cur_beta_theta[0])) = *(&(end_beta_theta[0]));
}

#endif
