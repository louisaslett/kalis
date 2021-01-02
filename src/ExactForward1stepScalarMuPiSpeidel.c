#include "ExactForward.h"
#include "Stencil.h"

#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_SCALAR
#define KALIS_1STEP
#define KALIS_SPEIDEL

#include "ExactForwardStencil.c"

#undef KALIS_MU
#undef KALIS_PI
#undef KALIS_1STEP
#undef KALIS_SPEIDEL
