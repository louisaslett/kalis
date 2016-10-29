#include <Rcpp.h>
using namespace Rcpp;

#include <stdlib.h>

#include "Cache.h"

// [[Rcpp::export]]
NumericMatrix ExactBackwardNaiveC_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
  NumericMatrix beta(N, N);
  double theta;
  char donor_hap, recipient_hap, recipient_hap_prev, H;
  int_fast32_t l=L-1;

  double *g, *gold;
  g    = (double*) malloc(sizeof(double)*N);
  gold = (double*) malloc(sizeof(double)*N);

  // Locus zero setup
  for(int_fast32_t recipient=0; recipient<N; ++recipient) {
    recipient_hap = seq_ind[recipient][l/8];
    gold[recipient] = 0.0;

    for(int_fast32_t donor=0; donor<N; ++donor) {
      donor_hap = seq_ind[donor][l/8];
      H = ((recipient_hap ^ donor_hap) >> (l%8)) & 1;
      theta = (H * mu[l]
                 + (1-H) * (1.0 - mu[l]));

      beta(donor, recipient) = 0.0;

      gold[recipient] += Pi(donor, recipient) * theta;
    }

    gold[recipient] = -log(gold[recipient]*rho[l-1]);
  }
  Rcout << "\n";

  while(l>t) {
    --l;
    for(int_fast32_t recipient=0; recipient<N; ++recipient) {
      recipient_hap_prev = seq_ind[recipient][(l+1)/8];
      recipient_hap = seq_ind[recipient][l/8];
      g[recipient] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        donor_hap = seq_ind[donor][(l+1)/8];
        H = ((recipient_hap_prev ^ donor_hap) >> ((l+1)%8)) & 1;
        theta = (H * mu[l+1]
                   + (1-H) * (1.0 - mu[l+1]));

        beta(donor, recipient) = (donor!=recipient) * (1.0 + (1.0 - rho[l]) * theta * exp(beta(donor, recipient) + gold[recipient]));

        donor_hap = seq_ind[donor][l/8];
        H = ((recipient_hap ^ donor_hap) >> (l%8)) & 1;
        theta = (H * mu[l]
                   + (1-H) * (1.0 - mu[l]));

        g[recipient] += Pi(donor, recipient) * theta * beta(donor, recipient);

        beta(donor, recipient) = log(beta(donor, recipient)) - gold[recipient];
      }

      gold[recipient] = -(log(g[recipient] * rho[l-1]) - gold[recipient]);
    }
  }

  free(g);
  free(gold);

  return(beta);
}
