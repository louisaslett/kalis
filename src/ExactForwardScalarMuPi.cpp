#include "Stencil.h"

#define EXACTFORWARDNOEXP ExactForwardNoExpAVX3_scmuPi
#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_SCALAR

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
