#include "Stencil.h"

#define EXACTFORWARDNOEXP ExactForward_speidel_scmu
#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_MATRIX
#define KALIS_SPEIDEL

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
#undef KALIS_SPEIDEL



#define EXACTFORWARDNOEXP ExactForward_scmu
#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_MATRIX

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
