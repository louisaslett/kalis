#include <Rcpp.h>
using namespace Rcpp;
#include <stdlib.h>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <utility>

#include "Cache.h"

// [[Rcpp::export]]
void ExactBackwardNaiveC_cpp(NumericMatrix beta,
                             NumericVector beta_g,
                             NumericVector beta_g2,
                             const int beta_from_rec,
                             const int beta_t,
                             const int t,
                             const int from_rec,
                             const int to_rec,
                             const int L,
                             const int N,
                             NumericMatrix Pi,
                             NumericVector mu,
                             NumericVector rho) {
  int_fast32_t l=beta_t;

  double *g, *gold, *goldold;
  g    = (double*) malloc(sizeof(double)*(to_rec-from_rec));
  gold = &(beta_g[0]);
  goldold = &(beta_g2[0]);

  // Locus L setup
  if(l>L-1) {
    l = L-1;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
      int32_t recipient_hap = (seq_locus[l][recipient/32] >> recipient%32) & 1;

      gold[recipient_beta] = 0.0;
      goldold[recipient_beta] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (seq_locus[l][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu[l]
                          + (1-H) * (1.0 - mu[l]));

        beta(donor, recipient_beta) = 1.0;

        gold[recipient_beta] += Pi(donor, recipient) * theta;
      }

      gold[recipient_beta] = -log(gold[recipient_beta]*rho[l-1]);
    }
  }

  while(l>t) {
    --l;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
      int32_t recipient_hap_prev = (seq_locus[l+1][recipient/32] >> recipient%32) & 1;
      int32_t recipient_hap = (seq_locus[l][recipient/32] >> recipient%32) & 1;

      g[recipient-from_rec] = 0.0;

      double gratio = exp(gold[recipient_beta]-goldold[recipient_beta]);
      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (seq_locus[l+1][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap_prev ^ donor_hap) & 1;
        double theta = (H * mu[l+1]
                          + (1-H) * (1.0 - mu[l+1]));

        beta(donor, recipient_beta) = (donor!=recipient) * (1.0 + (1.0 - rho[l]) * theta * beta(donor, recipient_beta) * gratio);

        donor_hap = (seq_locus[l][donor/32] >> donor%32) & 1;
        H = (recipient_hap ^ donor_hap) & 1;
        theta = (H * mu[l]
                   + (1-H) * (1.0 - mu[l]));

        g[recipient-from_rec] += Pi(donor, recipient) * theta * beta(donor, recipient_beta);
      }

      goldold[recipient_beta] = gold[recipient_beta];
      gold[recipient_beta] = -(log(g[recipient-from_rec] * rho[l-1]) - gold[recipient_beta]);
    }
  }

  free(g);
}

// [[Rcpp::export]]
void ParExactBackwardNaiveC_cpp(NumericMatrix beta,
                                NumericVector beta_g,
                                NumericVector beta_g2,
                                const int beta_from_rec,
                                const int beta_t,
                                const int t,
                                const int from_rec,
                                const int to_rec,
                                const int L,
                                const int N,
                                NumericMatrix Pi,
                                NumericVector mu,
                                NumericVector rho,
                                const int nthreads) {
  std::vector<std::thread> threads;

  if(nthreads < 2) {
    Rcout << "Only use parallel function with at least 2 threads";
  }

  // round(seq(0, nthreads) * double(to_rec-from_rec) / double(nthreads)) + from_rec;
  double spacing = double(to_rec-from_rec) / double(nthreads);

  for(int_fast32_t i=0; i<nthreads; ++i) {
    threads.push_back(std::thread(ExactBackwardNaiveC_cpp, beta, beta_g, beta_g2, beta_from_rec, beta_t, t, round(from_rec + i*spacing), round(from_rec + (i+1)*spacing), L, N, Pi, mu, rho));
    // Rcout << "From: " << round(from_rec + i*spacing) << ", To: " << round(from_rec + (i+1)*spacing) << "\n";
  }

  for(auto& th : threads) {
    th.join();
  }
}

// [[Rcpp::export]]
void ExactBackwardNoExpAVX3_cpp(NumericMatrix beta,
                                NumericVector beta_g,
                                NumericVector beta_g2,
                                const int beta_from_rec,
                                const int beta_t,
                                const int t,
                                const int from_rec,
                                const int to_rec,
                                const int L,
                                const int N,
                                NumericMatrix Pi,
                                NumericVector mu,
                                NumericVector rho) {
  int_fast32_t l=beta_t;

  double *__restrict__ g, *__restrict__ gold, *__restrict__ goldold;
  g       = (double*) malloc(sizeof(double)*(to_rec-from_rec));
  gold    = &(beta_g[0]);
  goldold = &(beta_g2[0]);

  // Locus L setup
  if(l>L-1) {
    l = L-1;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
      int32_t recipient_hap = (seq_locus[l][recipient/32] >> recipient%32) & 1;

      gold[recipient_beta] = 0.0;
      goldold[recipient_beta] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (seq_locus[l][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu[l]
                          + (1-H) * (1.0 - mu[l]));

        beta(donor, recipient_beta) = 1.0;

        gold[recipient_beta] += Pi(donor, recipient) * theta;
      }

      gold[recipient_beta] = -log(gold[recipient_beta]*rho[l-1]);
    }
  }

  const int_fast32_t reset_l = l;
  __m256d _one = _mm256_set1_pd(1.0);

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
      recipient_hap_prev -= (seq_locus[l+1][recipient/32] >> recipient%32) & 1;
      recipient_hap      -= (seq_locus[l][recipient/32] >> recipient%32) & 1;
      __m256i _recipient_hap_prev = _mm256_set1_epi32(recipient_hap_prev);
      __m256i _recipient_hap      = _mm256_set1_epi32(recipient_hap);

      g[recipient-from_rec] = 0.0; // For accumulating scalar part in ragged end ...
      __m256d _g = _mm256_set1_pd(0.0); // ... and for vector part.

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ betaRow = &(beta[N*recipient_beta]);
      double *__restrict__ PiRow = &(Pi[N*recipient]);
      // The diagonal should be zero ... we can fix that after the loop, but in order to avoid accumulating
      // incorrectly in g, we want to force the diagonal of Pi to zero
      PiRow[recipient] = 0.0;

      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      double muTmp1a = 2.0 * mu[l+1] - 1.0, muTmp2a = - mu[l+1] + 1.0;
      __m256d _muTmp1a = _mm256_broadcast_sd(&muTmp1a), _muTmp2a = _mm256_broadcast_sd(&muTmp2a);
      double muTmp1b = 2.0 * mu[l] - 1.0, muTmp2b = - mu[l] + 1.0;
      __m256d _muTmp1b = _mm256_broadcast_sd(&muTmp1b), _muTmp2b = _mm256_broadcast_sd(&muTmp2b);
      for(int_fast32_t donoroff=0; donoroff<N/(32*8); ++donoroff) {
        // Load next 256 donors and XOR with recipients
        __m256i _HA = _mm256_xor_si256(_recipient_hap_prev, _mm256_load_si256((__m256i*) &(seq_locus[l+1][donoroff*8])));
        uint32_t *HA = (uint32_t*) &_HA;
        __m256i _HB = _mm256_xor_si256(_recipient_hap,      _mm256_load_si256((__m256i*) &(seq_locus[l][donoroff*8])));
        uint32_t *HB = (uint32_t*) &_HB;

        const uint32_t mask = 16843009;
        for(int_fast32_t donor=0; donor<((32*8)/4)/4; ++donor) {
          // IACA_START
          double *betaNow1 = betaRow + donoroff*32*8 + donor*4*4;
          double *betaNow2 = betaRow + donoroff*32*8 + donor*4*4 + 4;
          double *betaNow3 = betaRow + donoroff*32*8 + donor*4*4 + 8;
          double *betaNow4 = betaRow + donoroff*32*8 + donor*4*4 + 12;

          __m256d _beta1  = _mm256_loadu_pd(betaNow1);
          __m256d _beta2  = _mm256_loadu_pd(betaNow2);
          __m256d _beta3  = _mm256_loadu_pd(betaNow3);
          __m256d _beta4  = _mm256_loadu_pd(betaNow4);

          __m256d _theta1 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4)/8]) >> (((donor*4)%8)*4), mask))));
          __m256d _theta2 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+1)/8]) >> (((donor*4+1)%8)*4), mask))));
          __m256d _theta3 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+2)/8]) >> (((donor*4+2)%8)*4), mask))));
          __m256d _theta4 = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HA[(donor*4+3)/8]) >> (((donor*4+3)%8)*4), mask))));

          _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
          _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
          _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
          _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1

          _beta1          = _mm256_mul_pd(_beta1, _theta1); // (theta*beta)
          _beta2          = _mm256_mul_pd(_beta2, _theta2); // (theta*beta)
          _beta3          = _mm256_mul_pd(_beta3, _theta3); // (theta*beta)
          _beta4          = _mm256_mul_pd(_beta4, _theta4); // (theta*beta)

          _beta1          = _mm256_fmadd_pd(_beta1, _gratioMulOmRho, _one); // (1.0 + {theta*beta} * {(1-rho)gratio})
          _beta2          = _mm256_fmadd_pd(_beta2, _gratioMulOmRho, _one); // (1.0 + {theta*beta} * {(1-rho)gratio})
          _beta3          = _mm256_fmadd_pd(_beta3, _gratioMulOmRho, _one); // (1.0 + {theta*beta} * {(1-rho)gratio})
          _beta4          = _mm256_fmadd_pd(_beta4, _gratioMulOmRho, _one); // (1.0 + {theta*beta} * {(1-rho)gratio})

          __m256d _pi1    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4);
          __m256d _pi2    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 4);
          __m256d _pi3    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 8);
          __m256d _pi4    = _mm256_loadu_pd(PiRow    + donoroff*32*8 + donor*4*4 + 12);

          _theta1         = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HB[(donor*4)/8]) >> (((donor*4)%8)*4), mask))));
          _theta2         = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HB[(donor*4+1)/8]) >> (((donor*4+1)%8)*4), mask))));
          _theta3         = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HB[(donor*4+2)/8]) >> (((donor*4+2)%8)*4), mask))));
          _theta4         = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((HB[(donor*4+3)/8]) >> (((donor*4+3)%8)*4), mask))));

          _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
          _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
          _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
          _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1

          _theta1         = _mm256_mul_pd(_beta1, _theta1); // (theta*beta)
          _theta2         = _mm256_mul_pd(_beta2, _theta2); // (theta*beta)
          _theta3         = _mm256_mul_pd(_beta3, _theta3); // (theta*beta)
          _theta4         = _mm256_mul_pd(_beta4, _theta4); // (theta*beta)

          _g              = _mm256_fmadd_pd(_pi1, _theta1, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_pi2, _theta2, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_pi3, _theta3, _g); // g += Pi * {theta*beta}
          _g              = _mm256_fmadd_pd(_pi4, _theta4, _g); // g += Pi * {theta*beta}

          _mm256_storeu_pd(betaNow1, _beta1);
          _mm256_storeu_pd(betaNow2, _beta2);
          _mm256_storeu_pd(betaNow3, _beta3);
          _mm256_storeu_pd(betaNow4, _beta4);
        }
        // IACA_END
      }
      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=0; donor<N%(32*8); ++donor) {
        int32_t donor_hap = (seq_locus[l+1][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
        int32_t H = (recipient_hap_prev ^ donor_hap) & 1;
        betaRow[(N/(32*8))*32*8+donor] = 1.0 + (H * muTmp1a + muTmp2a) * betaRow[(N/(32*8))*32*8+donor] * gratioMulOmRho;

        donor_hap = (seq_locus[l][(N/(32*8))*8 + donor/32] >> (donor%32)) & 1;
        H = (recipient_hap ^ donor_hap) & 1;
        g[recipient-from_rec] += PiRow[(N/(32*8))*32*8+donor] * (H * muTmp1b + muTmp2b) * betaRow[(N/(32*8))*32*8+donor];
      }
      // Fix the diagonal of beta ... note that because we set the diagonal of Pi
      // to zero we don't need to adjust g
      betaRow[recipient] = 0.0;

      // Accumulate row sum into g[recipient-from_rec]
      _g = _mm256_hadd_pd(_g, _g);
      g[recipient-from_rec] += ((double*)&_g)[0] + ((double*)&_g)[2];

      goldold[recipient_beta] = gold[recipient_beta];
      gold[recipient_beta] = -(log(g[recipient-from_rec] * rho[l-1]) - gold[recipient_beta]);
    }
  }

  free(g);
}

// [[Rcpp::export]]
void ParExactBackwardNoExpAVX3_cpp(NumericMatrix beta,
                                   NumericVector beta_g,
                                   NumericVector beta_g2,
                                   const int beta_from_rec,
                                   const int beta_t,
                                   const int t,
                                   const int from_rec,
                                   const int to_rec,
                                   const int L,
                                   const int N,
                                   NumericMatrix Pi,
                                   NumericVector mu,
                                   NumericVector rho,
                                   const int nthreads) {
  std::vector<std::thread> threads;

  if(nthreads < 2) {
    Rcout << "Only use parallel function with at least 2 threads";
  }

  // round(seq(0, nthreads) * double(to_rec-from_rec) / double(nthreads)) + from_rec;
  double spacing = double(to_rec-from_rec) / double(nthreads);

  for(int_fast32_t i=0; i<nthreads; ++i) {
    threads.push_back(std::thread(ExactBackwardNoExpAVX3_cpp, beta, beta_g, beta_g2, beta_from_rec, beta_t, t, round(from_rec + i*spacing), round(from_rec + (i+1)*spacing), L, N, Pi, mu, rho));
    // Rcout << "From: " << round(from_rec + i*spacing) << ", To: " << round(from_rec + (i+1)*spacing) << "\n";
  }

  for(auto& th : threads) {
    th.join();
  }
}

// [[Rcpp::export]]
void ExactBackwardNoExpAVX3single_cpp(NumericMatrix beta,
                                      NumericVector beta_g,
                                      NumericVector beta_g2,
                                      const int beta_from_rec,
                                      const int beta_t,
                                      const int t,
                                      const int from_rec,
                                      const int to_rec,
                                      const int L,
                                      const int N,
                                      NumericMatrix Pi,
                                      NumericVector mu,
                                      NumericVector rho) {
  int_fast32_t l=beta_t;

  double *__restrict__ g, *__restrict__ gold, *__restrict__ goldold;
  g       = (double*) malloc(sizeof(double)*(to_rec-from_rec));
  gold    = &(beta_g[0]);
  goldold = &(beta_g2[0]);
  double *Hcache_lp1, *Hcache_l, *Hcache_pointer_l, *Hcache_pointer_lp1;
  __m256d *_Hcache_lp1, *_Hcache_l;
  posix_memalign((void**) &Hcache_pointer_l, 32, sizeof(double)*N);
  posix_memalign((void**) &Hcache_pointer_lp1, 32, sizeof(double)*N);
  Hcache_lp1  = Hcache_pointer_lp1;
  Hcache_l    = Hcache_pointer_l;
  _Hcache_lp1 = (__m256d*) Hcache_lp1;
  _Hcache_l   = (__m256d*) Hcache_l;

  // Locus L setup
  if(l>L-1) {
    l = L-1;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;
      int32_t recipient_hap = (seq_locus[l][recipient/32] >> recipient%32) & 1;

      gold[recipient_beta] = 0.0;
      goldold[recipient_beta] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (seq_locus[l][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu[l]
                          + (1-H) * (1.0 - mu[l]));

        beta(donor, recipient_beta) = 1.0;

        gold[recipient_beta] += Pi(donor, recipient) * theta;
      }

      gold[recipient_beta] = -log(gold[recipient_beta]*rho[l-1]);
    }
  }

  const int_fast32_t reset_l = l;
  __m256d _one = _mm256_set1_pd(1.0);

  const uint32_t mask = 16843009;
  // Construct haplotype one ahead cache for this locus
  for(int_fast32_t indoff=0; indoff<N/32; indoff++) {
    _Hcache_l[indoff*8]     = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff], mask))));
    _Hcache_l[indoff*8 + 1] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 4, mask))));
    _Hcache_l[indoff*8 + 2] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 8, mask))));
    _Hcache_l[indoff*8 + 3] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 12, mask))));
    _Hcache_l[indoff*8 + 4] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 16, mask))));
    _Hcache_l[indoff*8 + 5] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 20, mask))));
    _Hcache_l[indoff*8 + 6] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 24, mask))));
    _Hcache_l[indoff*8 + 7] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 28, mask))));
  }
  for(int_fast32_t ind=(N/32)*32; ind<N; ind++) {
    Hcache_l[ind] = (double) ((seq_locus[l][ind/32] >> ind%32) & 1);
  }

  while(l>t) {
    --l;
    // Construct haplotype cache for this locus
    std::swap(Hcache_l, Hcache_lp1);
    std::swap(_Hcache_l, _Hcache_lp1);
    for(int_fast32_t indoff=0; indoff<N/32; indoff++) {
      _Hcache_l[indoff*8]     = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff], mask))));
      _Hcache_l[indoff*8 + 1] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 4, mask))));
      _Hcache_l[indoff*8 + 2] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 8, mask))));
      _Hcache_l[indoff*8 + 3] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 12, mask))));
      _Hcache_l[indoff*8 + 4] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 16, mask))));
      _Hcache_l[indoff*8 + 5] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 20, mask))));
      _Hcache_l[indoff*8 + 6] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 24, mask))));
      _Hcache_l[indoff*8 + 7] = _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32(seq_locus[l][indoff] >> 28, mask))));
    }
    for(int_fast32_t ind=(N/32)*32; ind<N; ind++) {
      Hcache_l[ind] = (double) ((seq_locus[l][ind/32] >> ind%32) & 1);
    }


    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_beta = recipient-beta_from_rec;

      // (1.0-rho[l]) * gratio .... for all recipients we consider to get us going
      double gratioMulOmRho = (1.0 - rho[l]) * exp(gold[recipient_beta] - goldold[recipient_beta]);
      __m256d _gratioMulOmRho = _mm256_set1_pd(gratioMulOmRho);

      // Load this recipient's bit into all 256-bits of an AVX register
      double recipient_hap_prev   = Hcache_lp1[recipient];
      double recipient_hap        = Hcache_l[recipient];
      __m256d _recipient_hap_prev = _mm256_set1_pd(Hcache_lp1[recipient]);
      __m256d _recipient_hap      = _mm256_set1_pd(Hcache_l[recipient]);

      g[recipient-from_rec] = 0.0; // For accumulating scalar part in ragged end ...
      __m256d _g = _mm256_set1_pd(0.0); // ... and for vector part.

      // TODO: for larger problems break this down into L1 cachable chunks of
      //       donors at a time
      double *__restrict__ betaRow = &(beta[N*recipient_beta]);
      double *__restrict__ PiRow = &(Pi[N*recipient]);
      // The diagonal should be zero ... we can fix that after the loop, but in order to avoid accumulating
      // incorrectly in g, we want to force the diagonal of Pi to zero
      PiRow[recipient] = 0.0;

      // Some temps to help computing (H * mu[l] + (1-H) * (1.0 - mu[l])) == H * (2*mu - 1) - mu + 1
      double muTmp1a = 2.0 * mu[l+1] - 1.0, muTmp2a = - mu[l+1] + 1.0;
      __m256d _muTmp1a = _mm256_broadcast_sd(&muTmp1a), _muTmp2a = _mm256_broadcast_sd(&muTmp2a);
      double muTmp1b = 2.0 * mu[l] - 1.0, muTmp2b = - mu[l] + 1.0;
      __m256d _muTmp1b = _mm256_broadcast_sd(&muTmp1b), _muTmp2b = _mm256_broadcast_sd(&muTmp2b);
      for(int_fast32_t donoroff=0; donoroff<N/(4*4); ++donoroff) {
        // IACA_START
        double *betaNow1 = betaRow + donoroff*4*4;
        double *betaNow2 = betaRow + donoroff*4*4 + 4;
        double *betaNow3 = betaRow + donoroff*4*4 + 8;
        double *betaNow4 = betaRow + donoroff*4*4 + 12;

        __m256d _beta1  = _mm256_loadu_pd(betaNow1);
        __m256d _beta2  = _mm256_loadu_pd(betaNow2);
        __m256d _beta3  = _mm256_loadu_pd(betaNow3);
        __m256d _beta4  = _mm256_loadu_pd(betaNow4);

        __m256d _theta1 = _mm256_xor_pd(_recipient_hap_prev, _Hcache_lp1[donoroff*4]);
        __m256d _theta2 = _mm256_xor_pd(_recipient_hap_prev, _Hcache_lp1[donoroff*4 + 1]);
        __m256d _theta3 = _mm256_xor_pd(_recipient_hap_prev, _Hcache_lp1[donoroff*4 + 2]);
        __m256d _theta4 = _mm256_xor_pd(_recipient_hap_prev, _Hcache_lp1[donoroff*4 + 3]);

        _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
        _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
        _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
        _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1

        _beta1          = _mm256_mul_pd(_beta1, _theta1); // (theta*beta)
        _beta2          = _mm256_mul_pd(_beta2, _theta2); // (theta*beta)
        _beta3          = _mm256_mul_pd(_beta3, _theta3); // (theta*beta)
        _beta4          = _mm256_mul_pd(_beta4, _theta4); // (theta*beta)

        _beta1          = _mm256_fmadd_pd(_beta1, _gratioMulOmRho, _one); // (1.0 + {theta*beta} * {(1-rho)gratio})
        _beta2          = _mm256_fmadd_pd(_beta2, _gratioMulOmRho, _one); // (1.0 + {theta*beta} * {(1-rho)gratio})
        _beta3          = _mm256_fmadd_pd(_beta3, _gratioMulOmRho, _one); // (1.0 + {theta*beta} * {(1-rho)gratio})
        _beta4          = _mm256_fmadd_pd(_beta4, _gratioMulOmRho, _one); // (1.0 + {theta*beta} * {(1-rho)gratio})

        __m256d _pi1    = _mm256_loadu_pd(PiRow + donoroff*4*4);
        __m256d _pi2    = _mm256_loadu_pd(PiRow + donoroff*4*4 + 4);
        __m256d _pi3    = _mm256_loadu_pd(PiRow + donoroff*4*4 + 8);
        __m256d _pi4    = _mm256_loadu_pd(PiRow + donoroff*4*4 + 12);

        _theta1         = _mm256_xor_pd(_recipient_hap, _Hcache_l[donoroff*4]);
        _theta2         = _mm256_xor_pd(_recipient_hap, _Hcache_l[donoroff*4 + 1]);
        _theta3         = _mm256_xor_pd(_recipient_hap, _Hcache_l[donoroff*4 + 2]);
        _theta4         = _mm256_xor_pd(_recipient_hap, _Hcache_l[donoroff*4 + 3]);

        _theta1         = _mm256_fmadd_pd(_theta1, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
        _theta2         = _mm256_fmadd_pd(_theta2, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
        _theta3         = _mm256_fmadd_pd(_theta3, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
        _theta4         = _mm256_fmadd_pd(_theta4, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1

        _theta1         = _mm256_mul_pd(_beta1, _theta1); // (theta*beta)
        _theta2         = _mm256_mul_pd(_beta2, _theta2); // (theta*beta)
        _theta3         = _mm256_mul_pd(_beta3, _theta3); // (theta*beta)
        _theta4         = _mm256_mul_pd(_beta4, _theta4); // (theta*beta)

        _g              = _mm256_fmadd_pd(_pi1, _theta1, _g); // g += Pi * {theta*beta}
        _g              = _mm256_fmadd_pd(_pi2, _theta2, _g); // g += Pi * {theta*beta}
        _g              = _mm256_fmadd_pd(_pi3, _theta3, _g); // g += Pi * {theta*beta}
        _g              = _mm256_fmadd_pd(_pi4, _theta4, _g); // g += Pi * {theta*beta}

        _mm256_storeu_pd(betaNow1, _beta1);
        _mm256_storeu_pd(betaNow2, _beta2);
        _mm256_storeu_pd(betaNow3, _beta3);
        _mm256_storeu_pd(betaNow4, _beta4);
      }
      // IACA_END

      // Tidy up any ragged end past a multiple of 256 ...
      for(int32_t donor=(N/(4*4))*4*4; donor<N; ++donor) {
        double donor_hap = Hcache_lp1[donor];
        double H = recipient_hap_prev + donor_hap - 2*recipient_hap_prev*donor_hap;
        betaRow[donor] = 1.0 + (H * muTmp1a + muTmp2a) * betaRow[donor] * gratioMulOmRho;

        donor_hap = Hcache_l[donor];
        H = recipient_hap + donor_hap - 2*recipient_hap*donor_hap;
        g[recipient-from_rec] += PiRow[donor] * (H * muTmp1b + muTmp2b) * betaRow[donor];
      }
      // Fix the diagonal of beta ... note that because we set the diagonal of Pi
      // to zero we don't need to adjust g
      betaRow[recipient] = 0.0;

      // Accumulate row sum into g[recipient-from_rec]
      _g = _mm256_hadd_pd(_g, _g);
      g[recipient-from_rec] += ((double*)&_g)[0] + ((double*)&_g)[2];

      goldold[recipient_beta] = gold[recipient_beta];
      gold[recipient_beta] = -(log(g[recipient-from_rec] * rho[l-1]) - gold[recipient_beta]);
    }
  }

  free(g);
  free(Hcache_pointer_l);
  free(Hcache_pointer_lp1);
}
