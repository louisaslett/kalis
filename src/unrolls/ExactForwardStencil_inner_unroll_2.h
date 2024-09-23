/* File generated by R/unroll-forward-inner-loop.R */

double *alphaNow0    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE);
double *alphaNow1    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (1*KALIS_DOUBLEVEC_SIZE);

KALIS_DOUBLE _alpha0 = KALIS_LOADU_DOUBLE(alphaNow0);
KALIS_DOUBLE _alpha1 = KALIS_LOADU_DOUBLE(alphaNow1);

KALIS_DOUBLE _pi0    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi1    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (1*KALIS_DOUBLEVEC_SIZE));

_pi0                 = KALIS_MUL_DOUBLE(_pi0, _rho);
_pi1                 = KALIS_MUL_DOUBLE(_pi1, _rho);

_alpha0              = KALIS_FMA_DOUBLE(_alpha0, _omRhoDivF, _pi0); // (Pi*rho + [(1-rho)/f] * alpha)
_alpha1              = KALIS_FMA_DOUBLE(_alpha1, _omRhoDivF, _pi1); // (Pi*rho + [(1-rho)/f] * alpha)

#if !defined(KALIS_1STEP)
KALIS_DOUBLE _theta0 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta1 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(1*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(1*KALIS_DOUBLEVEC_SIZE))%32));

_theta0              = KALIS_FMA_DOUBLE(_theta0, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta1              = KALIS_FMA_DOUBLE(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#else
KALIS_DOUBLE _theta0 = KALIS_SET_DOUBLE(1.0);
KALIS_DOUBLE _theta1 = KALIS_SET_DOUBLE(1.0);
#endif

_alpha0              = KALIS_MUL_DOUBLE(_theta0, _alpha0);
_alpha1              = KALIS_MUL_DOUBLE(_theta1, _alpha1);

_f                   = KALIS_ADD_DOUBLE(_f, _alpha0);
_f                   = KALIS_ADD_DOUBLE(_f, _alpha1);

KALIS_STOREU_DOUBLE(alphaNow0, _alpha0);
KALIS_STOREU_DOUBLE(alphaNow1, _alpha1);