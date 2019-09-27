#include <Rcpp.h>
using namespace Rcpp;

#include "Cache.h"
#include "ExactBackward.h"

// [[Rcpp::export]]
void ResetBackwardTable(List bck) {
  IntegerVector newl(1);
  newl[0] = 2147483647;
  bck["l"] = newl;
}

// [[Rcpp::export]]
void Backward_densePi_densemu_cpp(List bck,
                                  const int t,
                                  NumericMatrix Pi,
                                  NumericVector mu,
                                  NumericVector rho,
                                  const int nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix beta    = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g  = as<NumericVector>(bck["beta.g"]);
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
  if(nthreads>1) {
    ParExactBackwardNoExpAVX3_cpp(beta,
                                  beta_g,
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
  IntegerVector newl(1);
  newl[0] = t;
  bck["l"] = newl;
}

// [[Rcpp::export]]
void Backward_scalarPi_densemu_cpp(List bck,
                                   const int t,
                                   const double Pi,
                                   NumericVector mu,
                                   NumericVector rho,
                                   const int nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix beta    = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g  = as<NumericVector>(bck["beta.g"]);
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
  if(nthreads>1) {
    ParExactBackwardNoExpAVX3_scPi_cpp(beta,
                                       beta_g,
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
  IntegerVector newl(1);
  newl[0] = t;
  bck["l"] = newl;
}

// [[Rcpp::export]]
void Backward_densePi_scalarmu_cpp(List bck,
                                   const int t,
                                   NumericMatrix Pi,
                                   const double mu,
                                   NumericVector rho,
                                   const int nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix beta    = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g  = as<NumericVector>(bck["beta.g"]);
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
  if(nthreads>1) {
    ParExactBackwardNoExpAVX3_scmu_cpp(beta,
                                       beta_g,
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
    ExactBackwardNoExpAVX3_scmu_cpp(beta,
                                    beta_g,
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
  IntegerVector newl(1);
  newl[0] = t;
  bck["l"] = newl;
}

// [[Rcpp::export]]
void Backward_scalarPi_scalarmu_cpp(List bck,
                                    const int t,
                                    const double Pi,
                                    const double mu,
                                    NumericVector rho,
                                    const int nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix beta    = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g  = as<NumericVector>(bck["beta.g"]);
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
  if(nthreads>1) {
    ParExactBackwardNoExpAVX3_scmuPi_cpp(beta,
                                         beta_g,
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
    ExactBackwardNoExpAVX3_scmuPi_cpp(beta,
                                      beta_g,
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
  IntegerVector newl(1);
  newl[0] = t;
  bck["l"] = newl;
}
