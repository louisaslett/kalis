#include "R_ComputeStatus.h"

#include "StencilVec.h"

SEXP ComputeStatus() {
  SEXP res = PROTECT(Rf_allocVector(STRSXP, 1));

  const char *status = "" // "\n C version is " STRINGIFY_MACRO(__STDC_VERSION__) "\n"
#if defined(__x86_64__)
    "\nRunning in 64-bit mode using x86-64 architecture.\n"
#elif defined(__i686__)
    "\nRunning in 32-bit mode using i686 architecture.\n"
#elif defined(__i386__)
    "\nRunning in 32-bit mode using i386 architecture.\n"
#elif defined(__aarch64__)
    "\nRunning in 64-bit mode using ARM A64 architecture.\n"
#else
    "\nRunning on unknown architecture.\n"
#endif
    "Loops unrolled to depth " STRINGIFY_MACRO(KALIS_UNROLL) ".\n"
#if defined(KALIS_ISA_AVX512)
    "Currently using AVX-512, AVX2, SSE2 and BMI2 CPU instruction set extensions.\n"
#elif defined(KALIS_ISA_AVX2)
    "Currently using AVX2, AVX, SSE4.1, SSE2, FMA and BMI2 CPU instruction set extensions.\n"
#elif defined(KALIS_ISA_NEON)
    "Currently using ARM NEON and NEON FMA CPU instruction set extensions.\n"
#else
    "\nCurrently not using any special instruction sets (WARNING: poor performance likely).\n"
    "If this is unexpected (e.g. your CPU is Intel Haswell or newer architecture, or ARMv7+ with NEON support), then ensure that you are targeting the native architecture in compilation.  The easiest method is to add/change the following line in ~/.R/Makevars\n"
    "CFLAGS=-march=native -mtune=native -O3\n"
#endif
  ;
  SET_STRING_ELT(res, 0, Rf_mkChar(status));

  UNPROTECT(1);
  return(res);
}

SEXP VectorBitWidth() {
  return(Rf_ScalarInteger(KALIS_INTVEC_SIZE*32));
}
