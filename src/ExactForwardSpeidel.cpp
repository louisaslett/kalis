#include "Stencil.h"

#define EXACTFORWARDNOEXP ExactForward_speidel
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_MATRIX
#define KALIS_SPEIDEL

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
#undef KALIS_SPEIDEL
