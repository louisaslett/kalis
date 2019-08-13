#include <Rcpp.h>
using namespace Rcpp;

#include <string.h> // memcpy

#include "Cache.h"
#include "ExactForward.h"


// [[Rcpp::export]]
void ResetTable(List tbl) {
  as<IntegerVector>(tbl["l"])[0] = 0;
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
