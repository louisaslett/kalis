#include "Stencil.h"

#define EXACTFORWARDNOEXP ExactForward1stepNoExpAVX3_scmuPi
#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_SCALAR
#define KALIS_1STEP

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
