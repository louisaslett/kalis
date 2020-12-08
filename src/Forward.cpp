#include <Rcpp.h>
using namespace Rcpp;

#include "Cache.h"
#include "ExactForward.h"

// [[Rcpp::export]]
void ResetForwardTable(List fwd) {
  IntegerVector newl(1);
  newl[0] = 0;
  fwd["l"] = newl;
}

// [[Rcpp::export]]
void Forward_densePi_densemu_cpp(List fwd,
                                 const int t,
                                 NumericMatrix Pi,
                                 NumericVector mu,
                                 NumericVector rho,
                                 const bool use_speidel,
                                 IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix alpha    = as<NumericMatrix>(fwd["alpha"]);
  NumericVector alpha_f  = as<NumericVector>(fwd["alpha.f"]);
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

  int numthreads;
  if(nthreads.length() > 1) {
    numthreads = nthreads.length();
#if !defined(KALIS_PTHREAD_H)
    Rcout << "Thread affinity not supported on this platform, running on " << numthreads << " cores without setting affinity.\n";
#endif
  } else {
    numthreads = nthreads[0];
  }

  if(numthreads>1) {
    if(use_speidel)
      ParExactForward_speidel_cpp(alpha,
                                  alpha_f,
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
    else
      ParExactForward_cpp(alpha,
                          alpha_f,
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
    if(use_speidel)
      ExactForward_speidel_cpp(alpha,
                               alpha_f,
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
    else
      ExactForward_cpp(alpha,
                       alpha_f,
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
  fwd["l"] = newl;
}

// [[Rcpp::export]]
void Forward_scalarPi_densemu_cpp(List fwd,
                                  const int t,
                                  const double Pi,
                                  NumericVector mu,
                                  NumericVector rho,
                                  const bool use_speidel,
                                  IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix alpha    = as<NumericMatrix>(fwd["alpha"]);
  NumericVector alpha_f  = as<NumericVector>(fwd["alpha.f"]);
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

  int numthreads;
  if(nthreads.length() > 1) {
    numthreads = nthreads.length();
#if !defined(KALIS_PTHREAD_H)
    Rcout << "Thread affinity not supported on this platform, running on " << numthreads << " cores without setting affinity.\n";
#endif
  } else {
    numthreads = nthreads[0];
  }

  if(numthreads>1) {
    if(use_speidel)
      ParExactForward_speidel_scPi_cpp(alpha,
                                       alpha_f,
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
    else
      ParExactForward_scPi_cpp(alpha,
                               alpha_f,
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
    if(use_speidel)
      ExactForward_speidel_scPi_cpp(alpha,
                                    alpha_f,
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
    else
      ExactForward_scPi_cpp(alpha,
                            alpha_f,
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
  fwd["l"] = newl;
}

// [[Rcpp::export]]
void Forward_densePi_scalarmu_cpp(List fwd,
                                  const int t,
                                  NumericMatrix Pi,
                                  const double mu,
                                  NumericVector rho,
                                  const bool use_speidel,
                                  IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix alpha    = as<NumericMatrix>(fwd["alpha"]);
  NumericVector alpha_f  = as<NumericVector>(fwd["alpha.f"]);
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

  int numthreads;
  if(nthreads.length() > 1) {
    numthreads = nthreads.length();
#if !defined(KALIS_PTHREAD_H)
    Rcout << "Thread affinity not supported on this platform, running on " << numthreads << " cores without setting affinity.\n";
#endif
  } else {
    numthreads = nthreads[0];
  }

  if(numthreads>1) {
    if(use_speidel)
      ParExactForward_speidel_scmu_cpp(alpha,
                                       alpha_f,
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
    else
      ParExactForward_scmu_cpp(alpha,
                               alpha_f,
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
    if(use_speidel)
      ExactForward_speidel_scmu_cpp(alpha,
                                    alpha_f,
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
    else
      ExactForward_scmu_cpp(alpha,
                            alpha_f,
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
  fwd["l"] = newl;
}

// [[Rcpp::export]]
void Forward_scalarPi_scalarmu_cpp(List fwd,
                                   const int t,
                                   const double Pi,
                                   const double mu,
                                   NumericVector rho,
                                   const bool use_speidel,
                                   IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix alpha    = as<NumericMatrix>(fwd["alpha"]);
  NumericVector alpha_f  = as<NumericVector>(fwd["alpha.f"]);
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

  int numthreads;
  if(nthreads.length() > 1) {
    numthreads = nthreads.length();
#if !defined(KALIS_PTHREAD_H)
    Rcout << "Thread affinity not supported on this platform, running on " << numthreads << " cores without setting affinity.\n";
#endif
  } else {
    numthreads = nthreads[0];
  }

  if(numthreads>1) {
    if(use_speidel)
      ParExactForward_speidel_scmuPi_cpp(alpha,
                                         alpha_f,
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
    else
      ParExactForward_scmuPi_cpp(alpha,
                                 alpha_f,
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
    if(use_speidel)
      ExactForward_speidel_scmuPi_cpp(alpha,
                                      alpha_f,
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
    else
      ExactForward_scmuPi_cpp(alpha,
                              alpha_f,
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
  fwd["l"] = newl;
}

// [[Rcpp::export]]
void Forward1step_scalarPi_scalarmu_cpp(List fwd,
                                        const int t,
                                        const double Pi,
                                        const double mu,
                                        NumericVector rho,
                                        const bool use_speidel,
                                        IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix alpha    = as<NumericMatrix>(fwd["alpha"]);
  NumericVector alpha_f  = as<NumericVector>(fwd["alpha.f"]);
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

  int numthreads;
  if(nthreads.length() > 1) {
    numthreads = nthreads.length();
#if !defined(KALIS_PTHREAD_H)
    Rcout << "Thread affinity not supported on this platform, running on " << numthreads << " cores without setting affinity.\n";
#endif
  } else {
    numthreads = nthreads[0];
  }

  if(numthreads>1) {
    if(use_speidel)
      ParExactForward1step_speidel_scmuPi_cpp(alpha,
                                              alpha_f,
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
    else
      ParExactForward1step_scmuPi_cpp(alpha,
                                      alpha_f,
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
    if(use_speidel)
      ExactForward1step_speidel_scmuPi_cpp(alpha,
                                           alpha_f,
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
    else
      ExactForward1step_scmuPi_cpp(alpha,
                                   alpha_f,
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
  fwd["l"] = newl;
}
