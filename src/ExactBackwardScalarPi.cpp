#include "Stencil.h"

#define EXACTBACKWARDNOEXP ExactBackward_speidel_scPi
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_SCALAR
#define KALIS_SPEIDEL

#include "ExactBackwardStencil.cpp"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
#undef KALIS_SPEIDEL



#define EXACTBACKWARDNOEXP ExactBackward_scPi
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_SCALAR

#include "ExactBackwardStencil.cpp"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
