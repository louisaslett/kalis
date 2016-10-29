#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Exact_ComputeTable_naive_cpp(int l, NumericMatrix Pi, NumericVector mu, NumericVector rho) {
  return(Pi);
}
