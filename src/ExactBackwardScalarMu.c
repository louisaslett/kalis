#include "ExactBackward.h"
#include "Stencil.h"

#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_MATRIX

#include "ExactBackwardStencil.c"

#undef KALIS_MU
#undef KALIS_PI
