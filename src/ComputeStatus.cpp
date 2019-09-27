#include <Rcpp.h>
using namespace Rcpp;

#define KALIS_MU 1
#define KALIS_PI 1
#include "Stencil2.h"

// [[Rcpp::export]]
void ComputeStatus() {
#if defined(__x86_64__)
  Rcout << "\nRunning in 64-bit mode using x86-64 architecture.\n";
#elif defined(__i686__)
  Rcout << "\nRunning in 32-bit mode using i686 architecture.\n";
#elif defined(__i386__)
  Rcout << "\nRunning in 32-bit mode using i386 architecture.\n";
#else
  Rcout << "\nRunning on unknown architecture.\n";
#endif
  Rcout << "Loops unrolled to depth " << KALIS_UNROLL << ".\n";
#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)
  Rcout << "Currently using SSE2, SSE4.1, AVX, AVX2, FMA and BMI2 CPU instruction set extensions.\n\n";
#else
  Rcout << "\nCurrently not using any special instruction sets (WARNING: poor performance likely).\n";
  Rcout << "If this is unexpected (i.e. your CPU is Intel Haswell or newer architecture), then ensure that you are targeting the native architecture in compilation.  The easiest method is to add/change the following line in ~/.R/Makevars\n";
  Rcout << "CXX11FLAGS=-march=native -mtune=native -O3\n\n";
#endif
}
