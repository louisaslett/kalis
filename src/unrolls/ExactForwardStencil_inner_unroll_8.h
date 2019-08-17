/* File generated by R/unroll-forward-inner-loop.R */

double *alphaNow0    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE);
double *alphaNow1    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (1*KALIS_DOUBLEVEC_SIZE);
double *alphaNow2    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (2*KALIS_DOUBLEVEC_SIZE);
double *alphaNow3    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (3*KALIS_DOUBLEVEC_SIZE);
double *alphaNow4    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (4*KALIS_DOUBLEVEC_SIZE);
double *alphaNow5    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (5*KALIS_DOUBLEVEC_SIZE);
double *alphaNow6    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (6*KALIS_DOUBLEVEC_SIZE);
double *alphaNow7    = alphaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (7*KALIS_DOUBLEVEC_SIZE);

KALIS_DOUBLE _alpha0 = KALIS_LOADU_DOUBLE(alphaNow0);
KALIS_DOUBLE _alpha1 = KALIS_LOADU_DOUBLE(alphaNow1);
KALIS_DOUBLE _alpha2 = KALIS_LOADU_DOUBLE(alphaNow2);
KALIS_DOUBLE _alpha3 = KALIS_LOADU_DOUBLE(alphaNow3);
KALIS_DOUBLE _alpha4 = KALIS_LOADU_DOUBLE(alphaNow4);
KALIS_DOUBLE _alpha5 = KALIS_LOADU_DOUBLE(alphaNow5);
KALIS_DOUBLE _alpha6 = KALIS_LOADU_DOUBLE(alphaNow6);
KALIS_DOUBLE _alpha7 = KALIS_LOADU_DOUBLE(alphaNow7);

KALIS_DOUBLE _pi0    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi1    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (1*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi2    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (2*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi3    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (3*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi4    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (4*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi5    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (5*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi6    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (6*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi7    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (7*KALIS_DOUBLEVEC_SIZE));

_pi0                 = KALIS_MUL_DOUBLE(_pi0, _rho);
_pi1                 = KALIS_MUL_DOUBLE(_pi1, _rho);
_pi2                 = KALIS_MUL_DOUBLE(_pi2, _rho);
_pi3                 = KALIS_MUL_DOUBLE(_pi3, _rho);
_pi4                 = KALIS_MUL_DOUBLE(_pi4, _rho);
_pi5                 = KALIS_MUL_DOUBLE(_pi5, _rho);
_pi6                 = KALIS_MUL_DOUBLE(_pi6, _rho);
_pi7                 = KALIS_MUL_DOUBLE(_pi7, _rho);

_alpha0              = KALIS_FMA_DOUBLE(_alpha0, _omRhoDivF, _pi0); // (Pi*rho + [(1-rho)/f] * alpha)
_alpha1              = KALIS_FMA_DOUBLE(_alpha1, _omRhoDivF, _pi1); // (Pi*rho + [(1-rho)/f] * alpha)
_alpha2              = KALIS_FMA_DOUBLE(_alpha2, _omRhoDivF, _pi2); // (Pi*rho + [(1-rho)/f] * alpha)
_alpha3              = KALIS_FMA_DOUBLE(_alpha3, _omRhoDivF, _pi3); // (Pi*rho + [(1-rho)/f] * alpha)
_alpha4              = KALIS_FMA_DOUBLE(_alpha4, _omRhoDivF, _pi4); // (Pi*rho + [(1-rho)/f] * alpha)
_alpha5              = KALIS_FMA_DOUBLE(_alpha5, _omRhoDivF, _pi5); // (Pi*rho + [(1-rho)/f] * alpha)
_alpha6              = KALIS_FMA_DOUBLE(_alpha6, _omRhoDivF, _pi6); // (Pi*rho + [(1-rho)/f] * alpha)
_alpha7              = KALIS_FMA_DOUBLE(_alpha7, _omRhoDivF, _pi7); // (Pi*rho + [(1-rho)/f] * alpha)

#if !defined(KALIS_1STEP)
KALIS_DOUBLE _theta0 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta1 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(1*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(1*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta2 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(2*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(2*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta3 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(3*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(3*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta4 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(4*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(4*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta5 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(5*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(5*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta6 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(6*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(6*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta7 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(7*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(7*KALIS_DOUBLEVEC_SIZE))%32));

_theta0              = KALIS_FMA_DOUBLE(_theta0, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta1              = KALIS_FMA_DOUBLE(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta2              = KALIS_FMA_DOUBLE(_theta2, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta3              = KALIS_FMA_DOUBLE(_theta3, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta4              = KALIS_FMA_DOUBLE(_theta4, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta5              = KALIS_FMA_DOUBLE(_theta5, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta6              = KALIS_FMA_DOUBLE(_theta6, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta7              = KALIS_FMA_DOUBLE(_theta7, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#else
KALIS_DOUBLE _theta0 = KALIS_SET_DOUBLE(1.0);
KALIS_DOUBLE _theta1 = KALIS_SET_DOUBLE(1.0);
KALIS_DOUBLE _theta2 = KALIS_SET_DOUBLE(1.0);
KALIS_DOUBLE _theta3 = KALIS_SET_DOUBLE(1.0);
KALIS_DOUBLE _theta4 = KALIS_SET_DOUBLE(1.0);
KALIS_DOUBLE _theta5 = KALIS_SET_DOUBLE(1.0);
KALIS_DOUBLE _theta6 = KALIS_SET_DOUBLE(1.0);
KALIS_DOUBLE _theta7 = KALIS_SET_DOUBLE(1.0);
#endif

_alpha0              = KALIS_MUL_DOUBLE(_theta0, _alpha0);
_alpha1              = KALIS_MUL_DOUBLE(_theta1, _alpha1);
_alpha2              = KALIS_MUL_DOUBLE(_theta2, _alpha2);
_alpha3              = KALIS_MUL_DOUBLE(_theta3, _alpha3);
_alpha4              = KALIS_MUL_DOUBLE(_theta4, _alpha4);
_alpha5              = KALIS_MUL_DOUBLE(_theta5, _alpha5);
_alpha6              = KALIS_MUL_DOUBLE(_theta6, _alpha6);
_alpha7              = KALIS_MUL_DOUBLE(_theta7, _alpha7);

_f                   = KALIS_ADD_DOUBLE(_f, _alpha0);
_f                   = KALIS_ADD_DOUBLE(_f, _alpha1);
_f                   = KALIS_ADD_DOUBLE(_f, _alpha2);
_f                   = KALIS_ADD_DOUBLE(_f, _alpha3);
_f                   = KALIS_ADD_DOUBLE(_f, _alpha4);
_f                   = KALIS_ADD_DOUBLE(_f, _alpha5);
_f                   = KALIS_ADD_DOUBLE(_f, _alpha6);
_f                   = KALIS_ADD_DOUBLE(_f, _alpha7);

KALIS_STOREU_DOUBLE(alphaNow0, _alpha0);
KALIS_STOREU_DOUBLE(alphaNow1, _alpha1);
KALIS_STOREU_DOUBLE(alphaNow2, _alpha2);
KALIS_STOREU_DOUBLE(alphaNow3, _alpha3);
KALIS_STOREU_DOUBLE(alphaNow4, _alpha4);
KALIS_STOREU_DOUBLE(alphaNow5, _alpha5);
KALIS_STOREU_DOUBLE(alphaNow6, _alpha6);
KALIS_STOREU_DOUBLE(alphaNow7, _alpha7);
