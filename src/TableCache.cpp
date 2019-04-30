#include <Rcpp.h>
using namespace Rcpp;

#include <string.h> // memcpy

#include "Cache.h"
#include "ExactForward.h"

// [[Rcpp::export]]
void FillTableCache(List cache,
                    NumericMatrix Pi,
                    NumericVector mu,
                    NumericVector rho,
                    const int nthreads,
                    int from = 0,
                    int to = 0) {

  const int_fast32_t L = hap_size;
  if(from <= 0) {
    from = 1;
  }
  if(to <= 0 || to <= from || to >= L+1) {
    to = L;
  }
  const int_fast32_t chkpt_len = to-from;
  const int_fast32_t N = num_inds;
  const int_fast32_t cache_size = cache.length();
  const int_fast32_t alpha_size = as<NumericMatrix>(as<List>(cache[0])["alpha"]).length() * sizeof(double);
  const int_fast32_t alpha_f_size = as<NumericVector>(as<List>(cache[0])["alpha.f"]).length() * sizeof(double);
  const int_fast32_t alpha_f2_size = as<NumericVector>(as<List>(cache[0])["alpha.f2"]).length() * sizeof(double);

  double pos = 0.5;
  if(from > 1) {
    pos = 1.0;
  }
  as<NumericVector>(as<List>(cache[0])["l"])[0] = -1;
  for(int_fast32_t i = 0; i < cache_size; i++) {
    int_fast32_t t = floor((1.0-pos)*chkpt_len)+from-1;
    if(t == as<int>(as<List>(cache[i])["l"])-1) {
      as<NumericVector>(as<List>(cache[i])["l"])[0] = 0;
      break;
    }
    pos *= 0.5;

    Rcout << "Computing cache entry " << i+1 << " up to locus " << t+1
          << " for recipients " << as<int>(as<List>(cache[i])["from_recipient"])
          << " to " << as<int>(as<List>(cache[i])["to_recipient"])
          << " from ";
    if(as<int>(as<List>(cache[i])["l"]) < 1) {
      Rcout << "start\n";
    } else {
      Rcout << "locus " << as<int>(as<List>(cache[i])["l"]) << "\n";
    }

    if(nthreads>1) {
      ParExactForwardNoExpAVX3_cpp(as<NumericMatrix>(as<List>(cache[i])["alpha"]),
                                   as<NumericVector>(as<List>(cache[i])["alpha.f"]),
                                   as<NumericVector>(as<List>(cache[i])["alpha.f2"]),
                                   as<int>(as<List>(cache[i])["from_recipient"])-1,
                                   as<int>(as<List>(cache[i])["l"])-1,
                                   t,
                                   as<int>(as<List>(cache[i])["from_recipient"])-1,
                                   as<int>(as<List>(cache[i])["to_recipient"]),
                                   L,
                                   N,
                                   Pi,
                                   mu,
                                   rho,
                                   nthreads);
    } else {
      ExactForwardNoExpAVX3_cpp(as<NumericMatrix>(as<List>(cache[i])["alpha"]),
                                as<NumericVector>(as<List>(cache[i])["alpha.f"]),
                                as<NumericVector>(as<List>(cache[i])["alpha.f2"]),
                                as<int>(as<List>(cache[i])["from_recipient"])-1,
                                as<int>(as<List>(cache[i])["l"])-1,
                                t,
                                as<int>(as<List>(cache[i])["from_recipient"])-1,
                                as<int>(as<List>(cache[i])["to_recipient"]),
                                L,
                                N,
                                Pi,
                                mu,
                                rho);
    }
    as<NumericVector>(as<List>(cache[i])["l"])[0] = t+1;

    if(i < cache_size-1) { // Copy forward
      memcpy(&(as<NumericMatrix>(as<List>(cache[i+1])["alpha"])[0]),
             &(as<NumericMatrix>(as<List>(cache[i])["alpha"])[0]),
             alpha_size);
      memcpy(&(as<NumericVector>(as<List>(cache[i+1])["alpha.f"])[0]),
             &(as<NumericVector>(as<List>(cache[i])["alpha.f"])[0]),
             alpha_f_size);
      memcpy(&(as<NumericVector>(as<List>(cache[i+1])["alpha.f2"])[0]),
             &(as<NumericVector>(as<List>(cache[i])["alpha.f2"])[0]),
             alpha_f2_size);
      as<NumericVector>(as<List>(cache[i+1])["l"])[0] = as<NumericVector>(as<List>(cache[i])["l"])[0];
    }
  }
}

// [[Rcpp::export]]
void CopyForwardTable(List to, List from) {
  const int_fast32_t alpha_size = as<NumericMatrix>(to["alpha"]).length() * sizeof(double);
  const int_fast32_t alpha_f_size = as<NumericVector>(to["alpha.f"]).length() * sizeof(double);
  const int_fast32_t alpha_f2_size = as<NumericVector>(to["alpha.f2"]).length() * sizeof(double);

  memcpy(&(as<NumericMatrix>(to["alpha"])[0]),
         &(as<NumericMatrix>(from["alpha"])[0]),
         alpha_size);
  memcpy(&(as<NumericVector>(to["alpha.f"])[0]),
         &(as<NumericVector>(from["alpha.f"])[0]),
         alpha_f_size);
  memcpy(&(as<NumericVector>(to["alpha.f2"])[0]),
         &(as<NumericVector>(from["alpha.f2"])[0]),
         alpha_f2_size);
  as<IntegerVector>(to["l"])[0] = as<IntegerVector>(from["l"])[0];
}

// // [[Rcpp::export]]
// void ForwardUsingTableCache_cpp(List cache,
//                             const int t,
//                             NumericMatrix Pi,
//                             NumericVector mu,
//                             NumericVector rho,
//                             const int nthreads) {
//   const int L = hap_size;
//   const int N = num_inds;
//   NumericMatrix alpha    = as<NumericMatrix>(fwd["alpha"]);
//   NumericVector alpha_f  = as<NumericVector>(fwd["alpha.f"]);
//   NumericVector alpha_f2 = as<NumericVector>(fwd["alpha.f2"]);
//   int l        = as<int>(fwd["l"]);
//   int from_rec = as<int>(fwd["from_recipient"]);
//   int to_rec   = as<int>(fwd["to_recipient"]);
//
//   if(l>t) {
//     Rcout << "The forward table provided is for locus position " << l << " which is already past requested locus " << t << "\n";
//       return;
//   }
//   if(l==t) {
//     return;
//   }
//   if(nthreads>1) {
//     ParExactForwardNoExpAVX3_cpp(alpha,
//                                  alpha_f,
//                                  alpha_f2,
//                                  from_rec-1,
//                                  l-1,
//                                  t-1,
//                                  from_rec-1,
//                                  to_rec,
//                                  L,
//                                  N,
//                                  Pi,
//                                  mu,
//                                  rho,
//                                  nthreads);
//   } else {
//     ExactForwardNoExpAVX3_cpp(alpha,
//                               alpha_f,
//                               alpha_f2,
//                               from_rec-1,
//                               l-1,
//                               t-1,
//                               from_rec-1,
//                               to_rec,
//                               L,
//                               N,
//                               Pi,
//                               mu,
//                               rho);
//   }
//   as<NumericVector>(fwd["l"])[0] = t;
// }
