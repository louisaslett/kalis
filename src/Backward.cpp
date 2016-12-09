#include <Rcpp.h>
using namespace Rcpp;

#include "Cache.h"
#include "ExactBackward.h"

// [[Rcpp::export]]
void Backward(List bck,
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
  as<NumericVector>(bck["l"])[0] = t;
}
