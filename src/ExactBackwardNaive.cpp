#include <Rcpp.h>
using namespace Rcpp;
#include <stdlib.h>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <utility>
#include <functional>

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
