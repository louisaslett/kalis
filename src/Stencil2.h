#ifndef STENCIL2_H
#define STENCIL2_H

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

#define CPP_RAW_FN2(X) X ## _cpp_raw
#define CPP_RAW_FN(X) CPP_RAW_FN2(X)

#define CPP_FN2(X) X ## _cpp
#define CPP_FN(X) CPP_FN2(X)

#define PAR_CPP_FN2(X) Par ## X ## _cpp
#define PAR_CPP_FN(X) PAR_CPP_FN2(X)

#endif
