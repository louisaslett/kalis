#include "ExactBackward.h"
#include "Stencil.h"

#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_SCALAR

#include "ExactBackwardStencil.c"

#undef KALIS_MU
#undef KALIS_PI
