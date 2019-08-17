#ifndef STENCIL2_H
#define STENCIL2_H

#include "Stencil.h"

// Mu definitions

#ifndef KALIS_MU
#error "KALIS_MU not defined"
#endif
#if KALIS_MU == MU_SCALAR
#define MU_TYPE_C double
#define MU_TYPE_CPP const double
#define MU_ARG_CPP mu
#elif KALIS_MU == MU_VECTOR
#define MU_TYPE_C double *const __restrict__
#define MU_TYPE_CPP NumericVector
#define MU_ARG_CPP &(mu[0])
#else
#error "The type of KALIS_MU is not recognised"
#endif

// Pi definitions

#ifndef KALIS_PI
#error "KALIS_PI not defined"
#endif
#if KALIS_PI == PI_SCALAR
#define PI_TYPE_C double
#define PI_TYPE_CPP const double
#define PI_ARG_CPP Pi
#elif KALIS_PI == PI_MATRIX
#define PI_TYPE_C double *const __restrict__
#define PI_TYPE_CPP NumericMatrix
#define PI_ARG_CPP &(Pi[0])
#else
#error "The type of KALIS_PI is not recognised"
#endif

// Architecture definitions
// Only first set commented, apply to all

#define STRINGIFY_MACRO(x) STR(x)
#define STR(x) #x
#define EXPAND(x) x
#define CONCAT(n1, n2) STRINGIFY_MACRO(EXPAND(n1)EXPAND(n2))
#define KALIS_FORWARD_INNER_UNROLLED(N) STRINGIFY_MACRO(EXPAND(unrolls/ExactForwardStencil_inner_unroll_)EXPAND(N)EXPAND(.h))

#ifndef KALIS_UNROLL
#define KALIS_UNROLL 4
#endif

#if defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)

// How many fundamental types in a vector?
#define KALIS_INTVEC_SIZE 8
#define KALIS_DOUBLEVEC_SIZE 4

// Types
#define KALIS_DOUBLE __m256d
#define KALIS_INT32 __m256i

// Vector instructions
#define KALIS_SET_DOUBLE(X) _mm256_set1_pd(X)
#define KALIS_SET_INT32(X) _mm256_set1_epi32(X)
#define KALIS_LOAD_INT_VEC(X) _mm256_load_si256((__m256i*) &(X))
#define KALIS_XOR_INT(X, Y) _mm256_xor_si256(X, Y)
#define KALIS_LOADU_DOUBLE(X) _mm256_loadu_pd(X)
#define KALIS_STOREU_DOUBLE(X, Y) _mm256_storeu_pd(X, Y)
#define KALIS_ADD_DOUBLE(X, Y) _mm256_add_pd(X, Y)
#define KALIS_MUL_DOUBLE(X, Y) _mm256_mul_pd(X, Y)
#define KALIS_FMA_DOUBLE(X, Y, Z) _mm256_fmadd_pd(X, Y, Z)
#define KALIS_HSUM_DOUBLE(X) \
    (X) = _mm256_hadd_pd(X, _mm256_permute4x64_pd(X, 27)); \
    (X) = _mm256_hadd_pd(X, X);
#define KALIS_SPREADBITSTO_DOUBLE(X) _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((X), 16843009))));

#else

#define KALIS_INTVEC_SIZE 1
#define KALIS_DOUBLEVEC_SIZE 1

#define KALIS_DOUBLE double
#define KALIS_INT32 int32_t

#define KALIS_SET_DOUBLE(X) (X)
#define KALIS_SET_INT32(X) (X)
#define KALIS_LOAD_INT_VEC(X) (X)
#define KALIS_XOR_INT(X, Y) (X) ^ (Y)
#define KALIS_LOADU_DOUBLE(X) *(X)
#define KALIS_STOREU_DOUBLE(X, Y) *(X) = (Y)
#define KALIS_ADD_DOUBLE(X, Y) (X) + (Y)
#define KALIS_MUL_DOUBLE(X, Y) (X) * (Y)
#define KALIS_FMA_DOUBLE(X, Y, Z) (X) * (Y) + (Z)
#define KALIS_HSUM_DOUBLE(X) (X)
#define KALIS_SPREADBITSTO_DOUBLE(X) (double) ((X) & 1)

#endif


// Function definitions

#define CPP_RAW_FN2(X) X ## _cpp_raw
#define CPP_RAW_FN(X) CPP_RAW_FN2(X)

#define CPP_FN2(X) X ## _cpp
#define CPP_FN(X) CPP_FN2(X)

#define PAR_CPP_FN2(X) Par ## X ## _cpp
#define PAR_CPP_FN(X) PAR_CPP_FN2(X)

#endif
