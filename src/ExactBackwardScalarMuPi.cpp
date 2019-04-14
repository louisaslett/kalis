#include "Stencil.h"

#define EXACTBACKWARDNOEXP ExactBackwardNoExpAVX3_scmuPi
#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_SCALAR

#include "ExactBackwardStencil.cpp"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
