#include "Stencil.h"

#define EXACTFORWARDNOEXP ExactForward_speidel_scPi
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_SCALAR
#define KALIS_SPEIDEL

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
#undef KALIS_SPEIDEL



#define EXACTFORWARDNOEXP ExactForward_scPi
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_SCALAR

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
