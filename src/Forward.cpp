#include <Rcpp.h>
using namespace Rcpp;

#include "Cache.h"
#include "ExactForward.h"

// [[Rcpp::export]]
void Forward(List fwd,
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
  as<NumericVector>(fwd["l"])[0] = t;
}
