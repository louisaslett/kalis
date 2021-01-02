#ifndef STENCILVEC_H
#define STENCILVEC_H

#ifndef KALIS_UNROLL
#define KALIS_UNROLL 4
#endif

// Architecture definitions
// Only first set commented, apply to all

#if !defined(KALIS_NOASM) && defined(__SSE2__) && defined(__SSE4_1__) && defined(__AVX__) && defined(__AVX2__) && defined(__FMA__) && defined(__BMI2__)

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
#define KALIS_ANDNOT_INT(X, Y) _mm256_andnot_si256(Y, X)
#define KALIS_LOADU_DOUBLE(X) _mm256_loadu_pd(X)
#define KALIS_STOREU_DOUBLE(X, Y) _mm256_storeu_pd(X, Y)
#define KALIS_ADD_DOUBLE(X, Y) _mm256_add_pd(X, Y)
#define KALIS_MUL_DOUBLE(X, Y) _mm256_mul_pd(X, Y)
#define KALIS_FMA_DOUBLE(X, Y, Z) _mm256_fmadd_pd(X, Y, Z)
#define KALIS_HSUM_DOUBLE(X)                           \
(X) = _mm256_hadd_pd(X, _mm256_permute4x64_pd(X, 27)); \
(X) = _mm256_hadd_pd(X, X);
#define KALIS_SPREADBITSTO_DOUBLE(X) _mm256_cvtepi32_pd(_mm_cvtepi8_epi32(_mm_set_epi32(0, 0, 0, _pdep_u32((X), 16843009))));

#elif !defined(KALIS_NOASM) && defined(__ARM_NEON) && defined(__ARM_FEATURE_FMA)

#define KALIS_INTVEC_SIZE 4
#define KALIS_DOUBLEVEC_SIZE 2

#define KALIS_DOUBLE float64x2_t
#define KALIS_INT32 int32x4_t

#define KALIS_SET_DOUBLE(X) vdupq_n_f64((float64_t) X)
#define KALIS_SET_INT32(X) vdupq_n_s32(X)
#define KALIS_LOAD_INT_VEC(X) vld1q_s32((int32_t const *) &(X))
#define KALIS_XOR_INT(X, Y) veorq_s32(X, Y)
#define KALIS_ANDNOT_INT(X, Y) vandq_s32(X, vmvnq_s32(Y))
#define KALIS_LOADU_DOUBLE(X) vld1q_f64((float64_t const *) X)
#define KALIS_STOREU_DOUBLE(X, Y) vst1q_f64((float64_t*) X, Y)
#define KALIS_ADD_DOUBLE(X, Y) vaddq_f64(X, Y)
#define KALIS_MUL_DOUBLE(X, Y) vmulq_f64(X, Y)
#define KALIS_FMA_DOUBLE(X, Y, Z) vfmaq_f64(Z, X, Y)
#define KALIS_HSUM_DOUBLE(X)                             \
((double *) &X)[0] += ((double *) &X)[1];
#define KALIS_SPREADBITSTO_DOUBLE(X) ((float64x2_t) {(double) ((X)&1), (double) (((X)>>1)&1)});

#else

#define KALIS_INTVEC_SIZE 1
#define KALIS_DOUBLEVEC_SIZE 1

#define KALIS_DOUBLE double
#define KALIS_INT32 int32_t

#define KALIS_SET_DOUBLE(X) (X)
#define KALIS_SET_INT32(X) (X)
#define KALIS_LOAD_INT_VEC(X) (X)
#define KALIS_XOR_INT(X, Y) (X) ^ (Y)
#define KALIS_ANDNOT_INT(X, Y) (X) & ~(Y)
#define KALIS_LOADU_DOUBLE(X) *(X)
#define KALIS_STOREU_DOUBLE(X, Y) *(X) = (Y)
#define KALIS_ADD_DOUBLE(X, Y) (X) + (Y)
#define KALIS_MUL_DOUBLE(X, Y) (X) * (Y)
#define KALIS_FMA_DOUBLE(X, Y, Z) (X) * (Y) + (Z)
#define KALIS_HSUM_DOUBLE(X) (X = X)
#define KALIS_SPREADBITSTO_DOUBLE(X) (double) ((X) & 1)

#endif

// Utility macros
#define STRINGIFY_MACRO(x) STR(x)
#define STR(x) #x
#define EXPAND(x) x
#define CONCAT(n1, n2) STRINGIFY_MACRO(EXPAND(n1)EXPAND(n2))

#endif
