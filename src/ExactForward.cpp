#include "Stencil.h"

#define EXACTFORWARDNOEXP ExactForwardNoExpAVX3
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_MATRIX

#include "ExactForwardStencil.cpp"

#undef EXACTFORWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
