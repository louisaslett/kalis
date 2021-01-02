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
#define MU_TYPE_C double *const restrict
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
#define PI_TYPE_C double *const restrict
#define PI_TYPE_CPP NumericMatrix
#define PI_ARG_CPP &(Pi[0])
#else
#error "The type of KALIS_PI is not recognised"
#endif

// Architecture definitions

#include "StencilVec.h"

#define KALIS_FORWARD_INNER_UNROLLED(N) STRINGIFY_MACRO(EXPAND(unrolls/ExactForwardStencil_inner_unroll_)EXPAND(N)EXPAND(.h))
#define KALIS_BACKWARD_INNER_UNROLLED(N, V) STRINGIFY_MACRO(EXPAND(unrolls/ExactBackwardStencil_inner_unroll_)EXPAND(V)EXPAND(_)EXPAND(N)EXPAND(.h))


// Function definitions

#define FWD_CORE_ARGS2(OS,SP,MU,PI) forward ## OS ## SP ## MU ## PI ## _core_args
#define FWD_CORE_ARGS(OS,SP,MU,PI) FWD_CORE_ARGS2(OS,SP,MU,PI)

#define FWD_ARGS2(OS,SP,MU,PI) forward ## OS ## SP ## MU ## PI ## _args
#define FWD_ARGS(OS,SP,MU,PI) FWD_ARGS2(OS,SP,MU,PI)

#define FWD_RAW_FN2(OS,SP,MU,PI) ExactForward ## OS ## SP ## MU ## PI ## _raw
#define FWD_RAW_FN(OS,SP,MU,PI) FWD_RAW_FN2(OS,SP,MU,PI)

#define FWD_FN2(OS,SP,MU,PI) ExactForward ## OS ## SP ## MU ## PI
#define FWD_FN(OS,SP,MU,PI) FWD_FN2(OS,SP,MU,PI)

#define BCK_CORE_ARGS2(SP,MU,PI) backward ## SP ## MU ## PI ## _core_args
#define BCK_CORE_ARGS(SP,MU,PI) BCK_CORE_ARGS2(SP,MU,PI)

#define BCK_ARGS2(SP,MU,PI) backward ## SP ## MU ## PI ## _args
#define BCK_ARGS(SP,MU,PI) BCK_ARGS2(SP,MU,PI)

#define BCK_RAW_FN2(SP,MU,PI) ExactBackward ## SP ## MU ## PI ## _raw
#define BCK_RAW_FN(SP,MU,PI) BCK_RAW_FN2(SP,MU,PI)

#define BCK_FN2(SP,MU,PI) ExactBackward ## SP ## MU ## PI
#define BCK_FN(SP,MU,PI) BCK_FN2(SP,MU,PI)

#endif
