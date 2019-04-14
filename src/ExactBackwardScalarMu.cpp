#include "Stencil.h"

#define EXACTBACKWARDNOEXP ExactBackwardNoExpAVX3_scmu
#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_MATRIX

#include "ExactBackwardStencil.cpp"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
