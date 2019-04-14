#include "Stencil.h"

#define EXACTBACKWARDNOEXP ExactBackwardNoExpAVX3
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_MATRIX

#include "ExactBackwardStencil.cpp"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
