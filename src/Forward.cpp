#include <Rcpp.h>
using namespace Rcpp;

#include "Cache.h"
#include "ExactForward.h"

// [[Rcpp::export]]
void ResetForwardTable(List fwd) {
  as<NumericVector>(fwd["l"])[0] = 0;
}

// [[Rcpp::export]]
void Forward_densePi_cpp(List fwd,
                         const int t,
                         NumericMatrix Pi,
                         NumericVector mu,
                         NumericVector rho,
                         const int nthreads) {
  const int L = seq_size;
  const int N = num_inds;
  NumericMatrix alpha    = as<NumericMatrix>(fwd["alpha"]);
  NumericVector alpha_f  = as<NumericVector>(fwd["alpha.f"]);
  NumericVector alpha_f2 = as<NumericVector>(fwd["alpha.f2"]);
  int l        = as<int>(fwd["l"]);
  int from_rec = as<int>(fwd["from_recipient"]);
  int to_rec   = as<int>(fwd["to_recipient"]);

  if(l>t) {
    Rcout << "The forward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
    return;
  }
  if(l==t) {
    return;
  }
#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)
  if(nthreads>1) {
    ParExactForwardNoExpAVX3_cpp(alpha,
                                 alpha_f,
                                 alpha_f2,
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
    ExactForwardNoExpAVX3_cpp(alpha,
                              alpha_f,
                              alpha_f2,
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
    ParExactForwardNaiveC_cpp(alpha,
                                 alpha_f,
                                 alpha_f2,
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
    ExactForwardNaiveC_cpp(alpha,
                              alpha_f,
                              alpha_f2,
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

  as<NumericVector>(fwd["l"])[0] = t;
}

// [[Rcpp::export]]
void Forward_scalarPi_cpp(List fwd,
                          const int t,
                          const double Pi,
                          NumericVector mu,
                          NumericVector rho,
                          const int nthreads) {
  const int L = seq_size;
  const int N = num_inds;
  NumericMatrix alpha    = as<NumericMatrix>(fwd["alpha"]);
  NumericVector alpha_f  = as<NumericVector>(fwd["alpha.f"]);
  NumericVector alpha_f2 = as<NumericVector>(fwd["alpha.f2"]);
  int l        = as<int>(fwd["l"]);
  int from_rec = as<int>(fwd["from_recipient"]);
  int to_rec   = as<int>(fwd["to_recipient"]);

  if(l>t) {
    Rcout << "The forward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
    return;
  }
  if(l==t) {
    return;
  }
#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)
  if(nthreads>1) {
    ParExactForwardNoExpAVX3_scPi_cpp(alpha,
                                      alpha_f,
                                      alpha_f2,
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
    ExactForwardNoExpAVX3_scPi_cpp(alpha,
                                   alpha_f,
                                   alpha_f2,
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
  NumericMatrix Pimat(N, N, Pi);
  Pimat.fill_diag(0.0);

  if(nthreads>1) {
    ParExactForwardNaiveC_cpp(alpha,
                                 alpha_f,
                                 alpha_f2,
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
    ExactForwardNaiveC_cpp(alpha,
                              alpha_f,
                              alpha_f2,
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

  as<NumericVector>(fwd["l"])[0] = t;
}
