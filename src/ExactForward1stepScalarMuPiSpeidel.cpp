#include "Stencil.h"

#define EXACTFORWARDNOEXP ExactForward1step_speidel_scmuPi
#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_SCALAR
#define KALIS_1STEP
#define KALIS_SPEIDEL

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
#undef KALIS_1STEP
#undef KALIS_SPEIDEL
