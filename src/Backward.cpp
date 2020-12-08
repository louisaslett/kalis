#include <Rcpp.h>
using namespace Rcpp;

#include "Cache.h"
#include "ExactBackward.h"

// [[Rcpp::export]]
void ResetBackwardTable(List bck) {
  IntegerVector newl(1);
  newl[0] = 2147483647;
  bck["l"] = newl;
  LogicalVector newbt(1);
  newbt[0] = false;
  bck["beta.theta"] = newbt;
}

// [[Rcpp::export]]
void Backward_densePi_densemu_cpp(List bck,
                                  LogicalVector end_beta_theta,
                                  const int t,
                                  NumericMatrix Pi,
                                  NumericVector mu,
                                  NumericVector rho,
                                  const bool use_speidel,
                                  IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix beta           = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g         = as<NumericVector>(bck["beta.g"]);
  LogicalVector cur_beta_theta = as<LogicalVector>(bck["beta.theta"]);
  int l        = as<int>(bck["l"]);
  int from_rec = as<int>(bck["from_recipient"]);
  int to_rec   = as<int>(bck["to_recipient"]);

  if(l<t) {
    Rcout << "The backward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
    return;
  }
  if(l==t && !(!cur_beta_theta[0] && end_beta_theta[0])) {
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
      ParExactBackward_speidel_cpp(beta,
                                   beta_g,
                                   cur_beta_theta,
                                   end_beta_theta,
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
      ParExactBackward_cpp(beta,
                           beta_g,
                           cur_beta_theta,
                           end_beta_theta,
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
      ExactBackward_speidel_cpp(beta,
                                beta_g,
                                cur_beta_theta,
                                end_beta_theta,
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
      ExactBackward_cpp(beta,
                        beta_g,
                        cur_beta_theta,
                        end_beta_theta,
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
                                   LogicalVector end_beta_theta,
                                   const int t,
                                   const double Pi,
                                   NumericVector mu,
                                   NumericVector rho,
                                   const bool use_speidel,
                                   IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix beta           = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g         = as<NumericVector>(bck["beta.g"]);
  LogicalVector cur_beta_theta = as<LogicalVector>(bck["beta.theta"]);
  int l        = as<int>(bck["l"]);
  int from_rec = as<int>(bck["from_recipient"]);
  int to_rec   = as<int>(bck["to_recipient"]);

  if(l<t) {
    Rcout << "The backward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
    return;
  }
  if(l==t && !(!cur_beta_theta[0] && end_beta_theta[0])) {
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
      ParExactBackward_speidel_scPi_cpp(beta,
                                        beta_g,
                                        cur_beta_theta,
                                        end_beta_theta,
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
      ParExactBackward_scPi_cpp(beta,
                                beta_g,
                                cur_beta_theta,
                                end_beta_theta,
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
      ExactBackward_speidel_scPi_cpp(beta,
                                     beta_g,
                                     cur_beta_theta,
                                     end_beta_theta,
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
      ExactBackward_scPi_cpp(beta,
                             beta_g,
                             cur_beta_theta,
                             end_beta_theta,
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
                                   LogicalVector end_beta_theta,
                                   const int t,
                                   NumericMatrix Pi,
                                   const double mu,
                                   NumericVector rho,
                                   const bool use_speidel,
                                   IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix beta           = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g         = as<NumericVector>(bck["beta.g"]);
  LogicalVector cur_beta_theta = as<LogicalVector>(bck["beta.theta"]);
  int l        = as<int>(bck["l"]);
  int from_rec = as<int>(bck["from_recipient"]);
  int to_rec   = as<int>(bck["to_recipient"]);

  if(l<t) {
    Rcout << "The backward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
    return;
  }
  if(l==t && !(!cur_beta_theta[0] && end_beta_theta[0])) {
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
      ParExactBackward_speidel_scmu_cpp(beta,
                                        beta_g,
                                        cur_beta_theta,
                                        end_beta_theta,
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
      ParExactBackward_scmu_cpp(beta,
                                beta_g,
                                cur_beta_theta,
                                end_beta_theta,
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
      ExactBackward_speidel_scmu_cpp(beta,
                                     beta_g,
                                     cur_beta_theta,
                                     end_beta_theta,
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
      ExactBackward_scmu_cpp(beta,
                             beta_g,
                             cur_beta_theta,
                             end_beta_theta,
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
                                    LogicalVector end_beta_theta,
                                    const int t,
                                    const double Pi,
                                    const double mu,
                                    NumericVector rho,
                                    const bool use_speidel,
                                    IntegerVector nthreads) {
  const int L = hap_size;
  const int N = num_inds;
  NumericMatrix beta           = as<NumericMatrix>(bck["beta"]);
  NumericVector beta_g         = as<NumericVector>(bck["beta.g"]);
  LogicalVector cur_beta_theta = as<LogicalVector>(bck["beta.theta"]);
  int l        = as<int>(bck["l"]);
  int from_rec = as<int>(bck["from_recipient"]);
  int to_rec   = as<int>(bck["to_recipient"]);

  if(l<t) {
    Rcout << "The backward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
    return;
  }
  if(l==t && !(!cur_beta_theta[0] && end_beta_theta[0])) {
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
      ParExactBackward_speidel_scmuPi_cpp(beta,
                                          beta_g,
                                          cur_beta_theta,
                                          end_beta_theta,
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
      ParExactBackward_scmuPi_cpp(beta,
                                  beta_g,
                                  cur_beta_theta,
                                  end_beta_theta,
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
      ExactBackward_speidel_scmuPi_cpp(beta,
                                       beta_g,
                                       cur_beta_theta,
                                       end_beta_theta,
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
      ExactBackward_scmuPi_cpp(beta,
                               beta_g,
                               cur_beta_theta,
                               end_beta_theta,
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
