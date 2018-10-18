#include <Rcpp.h>
using namespace Rcpp;
#include <stdlib.h>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <functional>

#include "Cache.h"

// #include <iacaMarks.h>

// alpha: input (at locus alpha_t) and output (at locus t) forward matrix
// alpha_f: the vector of f's (at locus alpha_t) updated (to locus f)
// alpha_f2: the vector of f's (at locus alpha_t-1) updated (to locus f)
// alpha_from_rec: what recipient does row 1 of alpha, alpha_f and alpha_f2 correspond to?
// alpha_t: the zero indexed locus alpha is for, negative otherwise
// t: target zero indexed locus to compute up to
// from_rec/to_rec: zero indexed range of recipients to compute for ... [from,to)
// L: etc
// [[Rcpp::export]]
void ExactForwardNaiveC_cpp(NumericMatrix alpha,
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
                            NumericVector mu,
                            NumericVector rho) {
  int_fast32_t l=alpha_t;

  double *f, *fold, *foldold;
  f    = (double*) malloc(sizeof(double)*(to_rec-from_rec));
  fold = &(alpha_f[0]);
  foldold = &(alpha_f2[0]);

  // Locus zero setup
  if(l<0) {
    l = 0;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_alpha = recipient-alpha_from_rec;
      int32_t recipient_hap = (hap_locus[0][recipient/32] >> recipient%32) & 1;
      fold[recipient_alpha] = 0.0;
      foldold[recipient_alpha] = 0.0;

      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (hap_locus[0][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu[0]
                          + (1-H) * (1.0 - mu[0]));

        fold[recipient_alpha] += alpha(donor, recipient_alpha) = theta*Pi(donor, recipient);
      }

      fold[recipient_alpha] = -log(fold[recipient_alpha]*rho[0]);
    }
  }

  while(l<t) {
    ++l;
    for(int_fast32_t recipient=from_rec; recipient<to_rec; ++recipient) {
      int_fast32_t recipient_alpha = recipient-alpha_from_rec;
      int32_t recipient_hap = (hap_locus[l][recipient/32] >> recipient%32) & 1;

      f[recipient-from_rec] = 0.0;

      double fratio = exp(fold[recipient_alpha]-foldold[recipient_alpha]);
      for(int_fast32_t donor=0; donor<N; ++donor) {
        int32_t donor_hap = (hap_locus[l][donor/32] >> donor%32) & 1;
        int32_t H = (recipient_hap ^ donor_hap) & 1;
        double theta = (H * mu[l]
                          + (1-H) * (1.0 - mu[l]));

        f[recipient-from_rec] += alpha(donor, recipient_alpha) = (theta * Pi(donor, recipient)
                                                                    + theta * (1.0-rho[l-1]) * alpha(donor, recipient_alpha) * fratio);
      }

      foldold[recipient_alpha] = fold[recipient_alpha];
      fold[recipient_alpha] = -(log(f[recipient-from_rec] * rho[l]) - fold[recipient_alpha]);
    }
  }

  free(f);
}

// [[Rcpp::export]]
void ParExactForwardNaiveC_cpp(NumericMatrix alpha,
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
                               NumericVector mu,
                               NumericVector rho,
                               const int nthreads) {
  std::vector<std::thread> threads;

  if(nthreads < 2) {
    Rcout << "Only use parallel function with at least 2 threads";
  }

  // round(hap(0, nthreads) * double(to_rec-from_rec) / double(nthreads)) + from_rec;
  double spacing = double(to_rec-from_rec) / double(nthreads);

  for(int_fast32_t i=0; i<nthreads; ++i) {
    threads.push_back(std::thread(ExactForwardNaiveC_cpp, alpha, alpha_f, alpha_f2, alpha_from_rec, alpha_t, t, round(from_rec + i*spacing), round(from_rec + (i+1)*spacing), L, N, Pi, mu, rho));
    // Rcout << "From: " << round(from_rec + i*spacing) << ", To: " << round(from_rec + (i+1)*spacing) << "\n";
  }

  for(auto& th : threads) {
    th.join();
  }
}
