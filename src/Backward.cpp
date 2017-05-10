#include <Rcpp.h>
using namespace Rcpp;

#include "Cache.h"
#include "ExactBackward.h"

// [[Rcpp::export]]
void Backward_densePi_cpp(List bck,
                          const int t,
                          NumericMatrix Pi,
                          NumericVector mu,
                          NumericVector rho,
                          const int nthreads) {
  const int L = seq_size;
  const int N = num_inds;
  NumericMatrix beta    = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g  = as<NumericVector>(bck["beta.g"]);
  NumericVector beta_g2 = as<NumericVector>(bck["beta.g2"]);
  int l        = as<int>(bck["l"]);
  int from_rec = as<int>(bck["from_recipient"]);
  int to_rec   = as<int>(bck["to_recipient"]);

  if(l<t) {
    Rcout << "The backward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
    return;
  }
  if(l==t) {
    return;
  }
#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)
  if(nthreads>1) {
    ParExactBackwardNoExpAVX3_cpp(beta,
                                  beta_g,
                                  beta_g2,
                                  from_rec-1,
                                  l-1,
                                  t-1,
                                  from_rec-1,
                                  to_rec,
                                  L,
                                  N,
                                  Pi,
                                  mu,
                                  rho,
                                  nthreads);
  } else {
    ExactBackwardNoExpAVX3_cpp(beta,
                               beta_g,
                               beta_g2,
                               from_rec-1,
                               l-1,
                               t-1,
                               from_rec-1,
                               to_rec,
                               L,
                               N,
                               Pi,
                               mu,
                               rho);
  }
#else
  if(nthreads>1) {
    ParExactBackwardNaiveC_cpp(beta,
                                  beta_g,
                                  beta_g2,
                                  from_rec-1,
                                  l-1,
                                  t-1,
                                  from_rec-1,
                                  to_rec,
                                  L,
                                  N,
                                  Pi,
                                  mu,
                                  rho,
                                  nthreads);
  } else {
    ExactBackwardNaiveC_cpp(beta,
                               beta_g,
                               beta_g2,
                               from_rec-1,
                               l-1,
                               t-1,
                               from_rec-1,
                               to_rec,
                               L,
                               N,
                               Pi,
                               mu,
                               rho);
  }
#endif
  as<NumericVector>(bck["l"])[0] = t;
}

// [[Rcpp::export]]
void Backward_scalarPi_cpp(List bck,
                           const int t,
                           const double Pi,
                           NumericVector mu,
                           NumericVector rho,
                           const int nthreads) {
  const int L = seq_size;
  const int N = num_inds;
  NumericMatrix beta    = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g  = as<NumericVector>(bck["beta.g"]);
  NumericVector beta_g2 = as<NumericVector>(bck["beta.g2"]);
  int l        = as<int>(bck["l"]);
  int from_rec = as<int>(bck["from_recipient"]);
  int to_rec   = as<int>(bck["to_recipient"]);

  if(l<t) {
    Rcout << "The backward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
    return;
  }
  if(l==t) {
    return;
  }
#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)
  if(nthreads>1) {
    ParExactBackwardNoExpAVX3_scPi_cpp(beta,
                                       beta_g,
                                       beta_g2,
                                       from_rec-1,
                                       l-1,
                                       t-1,
                                       from_rec-1,
                                       to_rec,
                                       L,
                                       N,
                                       Pi,
                                       mu,
                                       rho,
                                       nthreads);
  } else {
    ExactBackwardNoExpAVX3_scPi_cpp(beta,
                                    beta_g,
                                    beta_g2,
                                    from_rec-1,
                                    l-1,
                                    t-1,
                                    from_rec-1,
                                    to_rec,
                                    L,
                                    N,
                                    Pi,
                                    mu,
                                    rho);
  }
#else
  NumericMatrix Pimat(N, N);
  Pimat.fill(Pi);
  Pimat.fill_diag(0.0);

  if(nthreads>1) {
    ParExactBackwardNaiveC_cpp(beta,
                                  beta_g,
                                  beta_g2,
                                  from_rec-1,
                                  l-1,
                                  t-1,
                                  from_rec-1,
                                  to_rec,
                                  L,
                                  N,
                                  Pimat,
                                  mu,
                                  rho,
                                  nthreads);
  } else {
    ExactBackwardNaiveC_cpp(beta,
                               beta_g,
                               beta_g2,
                               from_rec-1,
                               l-1,
                               t-1,
                               from_rec-1,
                               to_rec,
                               L,
                               N,
                               Pimat,
                               mu,
                               rho);
  }
#endif
  as<NumericVector>(bck["l"])[0] = t;
}
