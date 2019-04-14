#include "Stencil.h"

#define EXACTFORWARDNOEXP ExactForwardNoExpAVX3_scPi
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_SCALAR

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
