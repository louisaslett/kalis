#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void ComputeStatus() {
#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)
  Rcout << "\nCurrently using SSE2, SSE4.1, AVX, AVX2, FMA and BMI2 CPU instruction sets.\n\n";
#else
  Rcout << "\nCurrently not using any special instruction sets (WARNING: poor performance likely).\n";
  Rcout << "If this is unexpected (i.e. your CPU is Intel Haswell or newer architecture), then ensure that you are targeting the native architecture in compilation.  The easiest method is to add/change the following line in ~/.R/Makevars\n";
  Rcout << "CXX1XFLAGS=-march=native -mtune=native -O3\n\n";
#endif
}
