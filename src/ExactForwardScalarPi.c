#include "ExactForward.h"
#include "Stencil.h"

#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_SCALAR

#include "ExactForwardStencil.c"

#undef KALIS_MU
#undef KALIS_PI
