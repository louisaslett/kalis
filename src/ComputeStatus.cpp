#include <Rcpp.h>
using namespace Rcpp;

#define KALIS_MU 1
#define KALIS_PI 1
#include "Stencil2.h"

// [[Rcpp::export]]
std::string ComputeStatus() {
  std::string status;
#if defined(__x86_64__)
  status += "\nRunning in 64-bit mode using x86-64 architecture.\n";
#elif defined(__i686__)
  status += "\nRunning in 32-bit mode using i686 architecture.\n";
#elif defined(__i386__)
  status += "\nRunning in 32-bit mode using i386 architecture.\n";
#else
  status += "\nRunning on unknown architecture.\n";
#endif
  status += "Loops unrolled to depth " + std::to_string(KALIS_UNROLL) + ".\n";
#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)
  status += "Currently using SSE2, SSE4.1, AVX, AVX2, FMA and BMI2 CPU instruction set extensions.\n";
#else
  status += "\nCurrently not using any special instruction sets (WARNING: poor performance likely).\n";
  status += "If this is unexpected (i.e. your CPU is Intel Haswell or newer architecture), then ensure that you are targeting the native architecture in compilation.  The easiest method is to add/change the following line in ~/.R/Makevars\n";
  status += "CXX11FLAGS=-march=native -mtune=native -O3\n";
#endif
  return(status);
}
