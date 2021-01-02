#include "ExactBackward.h"
#include "Stencil.h"

#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_MATRIX

#include "ExactBackwardStencil.c"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
