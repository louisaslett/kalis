#include <Rcpp.h>
using namespace Rcpp;
#include <limits.h>

// [[Rcpp::export]]
void calc_tables_cpp_original_v1(NumericMatrix &cost_table,
                                 IntegerMatrix &index_table,
                                 NumericVector s,
                                 int maxcpt,
                                 int maxn) {
  for(int k = 0; k < maxcpt; k++) {
    for(int n = 0; n < maxn; n++) {
      double x = std::numeric_limits<double>::max();
      int xi;
      // Rcout << x << "\n";
      int y;
      for(int i = 0; i < n+1; i++) {
        y = cost_table(i, k+1) + s(i) + cost_table(n-i, k);
        // Rcout << cost_table(i,k+1) << " " << s(i) << " " << cost_table(n-i, k) << " = " << y << "\n";
        if(y < x) {
          xi = i;
          x = y;
        }
      }
      // Rcout << "Assign: " << x << " " << xi << "\n";
      index_table(n, k) = xi+1;
      cost_table(n+1, k+1) = x;
    }
  }
}


// [[Rcpp::export]]
void calc_tables_cpp_new_v2(NumericMatrix &cost_table,
                            IntegerMatrix &index_table,
                            NumericVector s,
                            int maxcpt,
                            int maxn) {
  for(int k = 0; k < maxcpt; k++) {
    for(int n = 0; n < maxn; n++) {
      double x = std::numeric_limits<double>::max();
      int xi;
      // Rcout << x << "\n";
      double y;
      for(int i = n; i > -1; i--) {
        y = cost_table(i, k+1) + s(i) + cost_table(n-i, k);
        // Rcout << cost_table(i,k+1) << " " << s(i) << " " << cost_table(n-i, k) << " = " << y << "\n";
        if(y <= x) {
          xi = i;
          x = y;
        } else {
          break;
        }
      }
      // Rcout << "Assign: " << x << " " << xi << "\n";
      index_table(n, k) = xi+1;
      cost_table(n+1, k+1) = x;
    }
  }
}
