#ifndef _PHASEGEOMF_F_H_
#define _PHASEGEOMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_GK_VELOCITY 
#define GUARDCOMPUTE_GK_VELOCITY 
// Prototype for Fortran procedure compute_gk_velocity ...
//
void FORTRAN_NAME( COMPUTE_GK_VELOCITY ,compute_gk_velocity )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(Z)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(larmor)
      ,CHFp_CONST_INT(include_drifts)
      ,CHFp_CONST_INT(include_par_streaming)
      ,CHFp_CONST_INT(include_gradb)
      ,CHFp_CONST_INT(mag_drifts_only)
      ,CHFp_CONST_INT(use_field_alignment)
      ,CHFp_CONST_FRA(E)
      ,CHFp_CONST_FRA(Bvec)
      ,CHFp_CONST_FRA(gradB)
      ,CHFp_CONST_FRA(curlb)
      ,CHFp_FRA(velocity) );

#define FORT_COMPUTE_GK_VELOCITY FORTRAN_NAME( inlineCOMPUTE_GK_VELOCITY, inlineCOMPUTE_GK_VELOCITY)
#define FORTNT_COMPUTE_GK_VELOCITY FORTRAN_NAME( COMPUTE_GK_VELOCITY, compute_gk_velocity)

inline void FORTRAN_NAME(inlineCOMPUTE_GK_VELOCITY, inlineCOMPUTE_GK_VELOCITY)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(Z)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(larmor)
      ,CHFp_CONST_INT(include_drifts)
      ,CHFp_CONST_INT(include_par_streaming)
      ,CHFp_CONST_INT(include_gradb)
      ,CHFp_CONST_INT(mag_drifts_only)
      ,CHFp_CONST_INT(use_field_alignment)
      ,CHFp_CONST_FRA(E)
      ,CHFp_CONST_FRA(Bvec)
      ,CHFp_CONST_FRA(gradB)
      ,CHFp_CONST_FRA(curlb)
      ,CHFp_FRA(velocity) )
{
 CH_TIMELEAF("FORT_COMPUTE_GK_VELOCITY");
 FORTRAN_NAME( COMPUTE_GK_VELOCITY ,compute_gk_velocity )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_REAL(Z)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_REAL(larmor)
      ,CHFt_CONST_INT(include_drifts)
      ,CHFt_CONST_INT(include_par_streaming)
      ,CHFt_CONST_INT(include_gradb)
      ,CHFt_CONST_INT(mag_drifts_only)
      ,CHFt_CONST_INT(use_field_alignment)
      ,CHFt_CONST_FRA(E)
      ,CHFt_CONST_FRA(Bvec)
      ,CHFt_CONST_FRA(gradB)
      ,CHFt_CONST_FRA(curlb)
      ,CHFt_FRA(velocity) );
}
#endif  // GUARDCOMPUTE_GK_VELOCITY 

#ifndef GUARDCOMPUTE_BFIELD_VELOCITY 
#define GUARDCOMPUTE_BFIELD_VELOCITY 
// Prototype for Fortran procedure compute_bfield_velocity ...
//
void FORTRAN_NAME( COMPUTE_BFIELD_VELOCITY ,compute_bfield_velocity )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(Bvec)
      ,CHFp_FRA(velocity) );

#define FORT_COMPUTE_BFIELD_VELOCITY FORTRAN_NAME( inlineCOMPUTE_BFIELD_VELOCITY, inlineCOMPUTE_BFIELD_VELOCITY)
#define FORTNT_COMPUTE_BFIELD_VELOCITY FORTRAN_NAME( COMPUTE_BFIELD_VELOCITY, compute_bfield_velocity)

inline void FORTRAN_NAME(inlineCOMPUTE_BFIELD_VELOCITY, inlineCOMPUTE_BFIELD_VELOCITY)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(Bvec)
      ,CHFp_FRA(velocity) )
{
 CH_TIMELEAF("FORT_COMPUTE_BFIELD_VELOCITY");
 FORTRAN_NAME( COMPUTE_BFIELD_VELOCITY ,compute_bfield_velocity )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(Bvec)
      ,CHFt_FRA(velocity) );
}
#endif  // GUARDCOMPUTE_BFIELD_VELOCITY 

#ifndef GUARDANNULUS_POLVEL_TEST 
#define GUARDANNULUS_POLVEL_TEST 
// Prototype for Fortran procedure annulus_polvel_test ...
//
void FORTRAN_NAME( ANNULUS_POLVEL_TEST ,annulus_polvel_test )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(rmin)
      ,CHFp_CONST_REAL(rbar)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_FRA(velocity)
      ,CHFp_INT(l_const_minorrad) );

#define FORT_ANNULUS_POLVEL_TEST FORTRAN_NAME( inlineANNULUS_POLVEL_TEST, inlineANNULUS_POLVEL_TEST)
#define FORTNT_ANNULUS_POLVEL_TEST FORTRAN_NAME( ANNULUS_POLVEL_TEST, annulus_polvel_test)

inline void FORTRAN_NAME(inlineANNULUS_POLVEL_TEST, inlineANNULUS_POLVEL_TEST)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(rmin)
      ,CHFp_CONST_REAL(rbar)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_FRA(velocity)
      ,CHFp_INT(l_const_minorrad) )
{
 CH_TIMELEAF("FORT_ANNULUS_POLVEL_TEST");
 FORTRAN_NAME( ANNULUS_POLVEL_TEST ,annulus_polvel_test )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_REAL(rmin)
      ,CHFt_CONST_REAL(rbar)
      ,CHFt_CONST_REAL(R0)
      ,CHFt_FRA(velocity)
      ,CHFt_INT(l_const_minorrad) );
}
#endif  // GUARDANNULUS_POLVEL_TEST 

#ifndef GUARDANNULUS_RADVEL_TEST 
#define GUARDANNULUS_RADVEL_TEST 
// Prototype for Fortran procedure annulus_radvel_test ...
//
void FORTRAN_NAME( ANNULUS_RADVEL_TEST ,annulus_radvel_test )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(rmin)
      ,CHFp_CONST_REAL(rbar)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_FRA(velocity)
      ,CHFp_INT(l_const_minorrad) );

#define FORT_ANNULUS_RADVEL_TEST FORTRAN_NAME( inlineANNULUS_RADVEL_TEST, inlineANNULUS_RADVEL_TEST)
#define FORTNT_ANNULUS_RADVEL_TEST FORTRAN_NAME( ANNULUS_RADVEL_TEST, annulus_radvel_test)

inline void FORTRAN_NAME(inlineANNULUS_RADVEL_TEST, inlineANNULUS_RADVEL_TEST)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(rmin)
      ,CHFp_CONST_REAL(rbar)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_FRA(velocity)
      ,CHFp_INT(l_const_minorrad) )
{
 CH_TIMELEAF("FORT_ANNULUS_RADVEL_TEST");
 FORTRAN_NAME( ANNULUS_RADVEL_TEST ,annulus_radvel_test )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_REAL(rmin)
      ,CHFt_CONST_REAL(rbar)
      ,CHFt_CONST_REAL(R0)
      ,CHFt_FRA(velocity)
      ,CHFt_INT(l_const_minorrad) );
}
#endif  // GUARDANNULUS_RADVEL_TEST 

#ifndef GUARDANNULUS_RADPOLVEL_TEST 
#define GUARDANNULUS_RADPOLVEL_TEST 
// Prototype for Fortran procedure annulus_radpolvel_test ...
//
void FORTRAN_NAME( ANNULUS_RADPOLVEL_TEST ,annulus_radpolvel_test )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(rmin)
      ,CHFp_CONST_REAL(rbar)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_FRA(velocity)
      ,CHFp_INT(l_const_minorrad) );

#define FORT_ANNULUS_RADPOLVEL_TEST FORTRAN_NAME( inlineANNULUS_RADPOLVEL_TEST, inlineANNULUS_RADPOLVEL_TEST)
#define FORTNT_ANNULUS_RADPOLVEL_TEST FORTRAN_NAME( ANNULUS_RADPOLVEL_TEST, annulus_radpolvel_test)

inline void FORTRAN_NAME(inlineANNULUS_RADPOLVEL_TEST, inlineANNULUS_RADPOLVEL_TEST)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(rmin)
      ,CHFp_CONST_REAL(rbar)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_FRA(velocity)
      ,CHFp_INT(l_const_minorrad) )
{
 CH_TIMELEAF("FORT_ANNULUS_RADPOLVEL_TEST");
 FORTRAN_NAME( ANNULUS_RADPOLVEL_TEST ,annulus_radpolvel_test )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_REAL(rmin)
      ,CHFt_CONST_REAL(rbar)
      ,CHFt_CONST_REAL(R0)
      ,CHFt_FRA(velocity)
      ,CHFt_INT(l_const_minorrad) );
}
#endif  // GUARDANNULUS_RADPOLVEL_TEST 

#ifndef GUARDMAJOR_RADIAL_VEL_TEST 
#define GUARDMAJOR_RADIAL_VEL_TEST 
// Prototype for Fortran procedure major_radial_vel_test ...
//
void FORTRAN_NAME( MAJOR_RADIAL_VEL_TEST ,major_radial_vel_test )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_FRA(velocity) );

#define FORT_MAJOR_RADIAL_VEL_TEST FORTRAN_NAME( inlineMAJOR_RADIAL_VEL_TEST, inlineMAJOR_RADIAL_VEL_TEST)
#define FORTNT_MAJOR_RADIAL_VEL_TEST FORTRAN_NAME( MAJOR_RADIAL_VEL_TEST, major_radial_vel_test)

inline void FORTRAN_NAME(inlineMAJOR_RADIAL_VEL_TEST, inlineMAJOR_RADIAL_VEL_TEST)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_FRA(velocity) )
{
 CH_TIMELEAF("FORT_MAJOR_RADIAL_VEL_TEST");
 FORTRAN_NAME( MAJOR_RADIAL_VEL_TEST ,major_radial_vel_test )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_FRA(velocity) );
}
#endif  // GUARDMAJOR_RADIAL_VEL_TEST 

#ifndef GUARDFREE_STREAM_VEL_TEST 
#define GUARDFREE_STREAM_VEL_TEST 
// Prototype for Fortran procedure free_stream_vel_test ...
//
void FORTRAN_NAME( FREE_STREAM_VEL_TEST ,free_stream_vel_test )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(component)
      ,CHFp_CONST_FRA1(major_radius)
      ,CHFp_CONST_INT(axisymmetric)
      ,CHFp_FRA(velocity) );

#define FORT_FREE_STREAM_VEL_TEST FORTRAN_NAME( inlineFREE_STREAM_VEL_TEST, inlineFREE_STREAM_VEL_TEST)
#define FORTNT_FREE_STREAM_VEL_TEST FORTRAN_NAME( FREE_STREAM_VEL_TEST, free_stream_vel_test)

inline void FORTRAN_NAME(inlineFREE_STREAM_VEL_TEST, inlineFREE_STREAM_VEL_TEST)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(component)
      ,CHFp_CONST_FRA1(major_radius)
      ,CHFp_CONST_INT(axisymmetric)
      ,CHFp_FRA(velocity) )
{
 CH_TIMELEAF("FORT_FREE_STREAM_VEL_TEST");
 FORTRAN_NAME( FREE_STREAM_VEL_TEST ,free_stream_vel_test )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(component)
      ,CHFt_CONST_FRA1(major_radius)
      ,CHFt_CONST_INT(axisymmetric)
      ,CHFt_FRA(velocity) );
}
#endif  // GUARDFREE_STREAM_VEL_TEST 

#ifndef GUARDCOMPUTE_BSTAR 
#define GUARDCOMPUTE_BSTAR 
// Prototype for Fortran procedure compute_bstar ...
//
void FORTRAN_NAME( COMPUTE_BSTAR ,compute_bstar )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(prefactor)
      ,CHFp_CONST_FRA(B)
      ,CHFp_CONST_FRA1(B_magnitude)
      ,CHFp_CONST_FRA(curlb)
      ,CHFp_CONST_FRA1(bdotcurlb)
      ,CHFp_FRA(BStar)
      ,CHFp_FRA1(BStarPar) );

#define FORT_COMPUTE_BSTAR FORTRAN_NAME( inlineCOMPUTE_BSTAR, inlineCOMPUTE_BSTAR)
#define FORTNT_COMPUTE_BSTAR FORTRAN_NAME( COMPUTE_BSTAR, compute_bstar)

inline void FORTRAN_NAME(inlineCOMPUTE_BSTAR, inlineCOMPUTE_BSTAR)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(prefactor)
      ,CHFp_CONST_FRA(B)
      ,CHFp_CONST_FRA1(B_magnitude)
      ,CHFp_CONST_FRA(curlb)
      ,CHFp_CONST_FRA1(bdotcurlb)
      ,CHFp_FRA(BStar)
      ,CHFp_FRA1(BStarPar) )
{
 CH_TIMELEAF("FORT_COMPUTE_BSTAR");
 FORTRAN_NAME( COMPUTE_BSTAR ,compute_bstar )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_REAL(prefactor)
      ,CHFt_CONST_FRA(B)
      ,CHFt_CONST_FRA1(B_magnitude)
      ,CHFt_CONST_FRA(curlb)
      ,CHFt_CONST_FRA1(bdotcurlb)
      ,CHFt_FRA(BStar)
      ,CHFt_FRA1(BStarPar) );
}
#endif  // GUARDCOMPUTE_BSTAR 

#ifndef GUARDPSCSPOINTDOTPRODCFG 
#define GUARDPSCSPOINTDOTPRODCFG 
// Prototype for Fortran procedure pscspointdotprodcfg ...
//
void FORTRAN_NAME( PSCSPOINTDOTPRODCFG ,pscspointdotprodcfg )(
      CHFp_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_INT(astart)
      ,CHFp_CONST_FRA(b)
      ,CHFp_CONST_INT(bstart)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_BOX(gridbox) );

#define FORT_PSCSPOINTDOTPRODCFG FORTRAN_NAME( inlinePSCSPOINTDOTPRODCFG, inlinePSCSPOINTDOTPRODCFG)
#define FORTNT_PSCSPOINTDOTPRODCFG FORTRAN_NAME( PSCSPOINTDOTPRODCFG, pscspointdotprodcfg)

inline void FORTRAN_NAME(inlinePSCSPOINTDOTPRODCFG, inlinePSCSPOINTDOTPRODCFG)(
      CHFp_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_INT(astart)
      ,CHFp_CONST_FRA(b)
      ,CHFp_CONST_INT(bstart)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_BOX(gridbox) )
{
 CH_TIMELEAF("FORT_PSCSPOINTDOTPRODCFG");
 FORTRAN_NAME( PSCSPOINTDOTPRODCFG ,pscspointdotprodcfg )(
      CHFt_FRA1(c)
      ,CHFt_CONST_FRA(a)
      ,CHFt_CONST_INT(astart)
      ,CHFt_CONST_FRA(b)
      ,CHFt_CONST_INT(bstart)
      ,CHFt_CONST_INT(ncomp)
      ,CHFt_BOX(gridbox) );
}
#endif  // GUARDPSCSPOINTDOTPRODCFG 

#ifndef GUARDPSCSPOINTDOTPRODVEL 
#define GUARDPSCSPOINTDOTPRODVEL 
// Prototype for Fortran procedure pscspointdotprodvel ...
//
void FORTRAN_NAME( PSCSPOINTDOTPRODVEL ,pscspointdotprodvel )(
      CHFp_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_INT(astart)
      ,CHFp_CONST_FRA(b)
      ,CHFp_CONST_INT(bstart)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_BOX(gridbox) );

#define FORT_PSCSPOINTDOTPRODVEL FORTRAN_NAME( inlinePSCSPOINTDOTPRODVEL, inlinePSCSPOINTDOTPRODVEL)
#define FORTNT_PSCSPOINTDOTPRODVEL FORTRAN_NAME( PSCSPOINTDOTPRODVEL, pscspointdotprodvel)

inline void FORTRAN_NAME(inlinePSCSPOINTDOTPRODVEL, inlinePSCSPOINTDOTPRODVEL)(
      CHFp_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_INT(astart)
      ,CHFp_CONST_FRA(b)
      ,CHFp_CONST_INT(bstart)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_BOX(gridbox) )
{
 CH_TIMELEAF("FORT_PSCSPOINTDOTPRODVEL");
 FORTRAN_NAME( PSCSPOINTDOTPRODVEL ,pscspointdotprodvel )(
      CHFt_FRA1(c)
      ,CHFt_CONST_FRA(a)
      ,CHFt_CONST_INT(astart)
      ,CHFt_CONST_FRA(b)
      ,CHFt_CONST_INT(bstart)
      ,CHFt_CONST_INT(ncomp)
      ,CHFt_BOX(gridbox) );
}
#endif  // GUARDPSCSPOINTDOTPRODVEL 

#ifndef GUARDPSCSPOINTWISEREDUCEDDOTPRODCFG 
#define GUARDPSCSPOINTWISEREDUCEDDOTPRODCFG 
// Prototype for Fortran procedure pscspointwisereduceddotprodcfg ...
//
void FORTRAN_NAME( PSCSPOINTWISEREDUCEDDOTPRODCFG ,pscspointwisereduceddotprodcfg )(
      CHFp_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_INT(astart)
      ,CHFp_CONST_FRA(b)
      ,CHFp_CONST_INT(bstart)
      ,CHFp_CONST_INT(numcomp)
      ,CHFp_BOX(gridbox) );

#define FORT_PSCSPOINTWISEREDUCEDDOTPRODCFG FORTRAN_NAME( inlinePSCSPOINTWISEREDUCEDDOTPRODCFG, inlinePSCSPOINTWISEREDUCEDDOTPRODCFG)
#define FORTNT_PSCSPOINTWISEREDUCEDDOTPRODCFG FORTRAN_NAME( PSCSPOINTWISEREDUCEDDOTPRODCFG, pscspointwisereduceddotprodcfg)

inline void FORTRAN_NAME(inlinePSCSPOINTWISEREDUCEDDOTPRODCFG, inlinePSCSPOINTWISEREDUCEDDOTPRODCFG)(
      CHFp_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_INT(astart)
      ,CHFp_CONST_FRA(b)
      ,CHFp_CONST_INT(bstart)
      ,CHFp_CONST_INT(numcomp)
      ,CHFp_BOX(gridbox) )
{
 CH_TIMELEAF("FORT_PSCSPOINTWISEREDUCEDDOTPRODCFG");
 FORTRAN_NAME( PSCSPOINTWISEREDUCEDDOTPRODCFG ,pscspointwisereduceddotprodcfg )(
      CHFt_FRA1(c)
      ,CHFt_CONST_FRA(a)
      ,CHFt_CONST_INT(astart)
      ,CHFt_CONST_FRA(b)
      ,CHFt_CONST_INT(bstart)
      ,CHFt_CONST_INT(numcomp)
      ,CHFt_BOX(gridbox) );
}
#endif  // GUARDPSCSPOINTWISEREDUCEDDOTPRODCFG 

#ifndef GUARDPSCSPOINTWISEREDUCEDDOTPRODVEL 
#define GUARDPSCSPOINTWISEREDUCEDDOTPRODVEL 
// Prototype for Fortran procedure pscspointwisereduceddotprodvel ...
//
void FORTRAN_NAME( PSCSPOINTWISEREDUCEDDOTPRODVEL ,pscspointwisereduceddotprodvel )(
      CHFp_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_INT(astart)
      ,CHFp_CONST_FRA(b)
      ,CHFp_CONST_INT(bstart)
      ,CHFp_CONST_INT(numcomp)
      ,CHFp_BOX(gridbox) );

#define FORT_PSCSPOINTWISEREDUCEDDOTPRODVEL FORTRAN_NAME( inlinePSCSPOINTWISEREDUCEDDOTPRODVEL, inlinePSCSPOINTWISEREDUCEDDOTPRODVEL)
#define FORTNT_PSCSPOINTWISEREDUCEDDOTPRODVEL FORTRAN_NAME( PSCSPOINTWISEREDUCEDDOTPRODVEL, pscspointwisereduceddotprodvel)

inline void FORTRAN_NAME(inlinePSCSPOINTWISEREDUCEDDOTPRODVEL, inlinePSCSPOINTWISEREDUCEDDOTPRODVEL)(
      CHFp_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_INT(astart)
      ,CHFp_CONST_FRA(b)
      ,CHFp_CONST_INT(bstart)
      ,CHFp_CONST_INT(numcomp)
      ,CHFp_BOX(gridbox) )
{
 CH_TIMELEAF("FORT_PSCSPOINTWISEREDUCEDDOTPRODVEL");
 FORTRAN_NAME( PSCSPOINTWISEREDUCEDDOTPRODVEL ,pscspointwisereduceddotprodvel )(
      CHFt_FRA1(c)
      ,CHFt_CONST_FRA(a)
      ,CHFt_CONST_INT(astart)
      ,CHFt_CONST_FRA(b)
      ,CHFt_CONST_INT(bstart)
      ,CHFt_CONST_INT(numcomp)
      ,CHFt_BOX(gridbox) );
}
#endif  // GUARDPSCSPOINTWISEREDUCEDDOTPRODVEL 

#ifndef GUARDPSCSTANGRADFACE 
#define GUARDPSCSTANGRADFACE 
// Prototype for Fortran procedure pscstangradface ...
//
void FORTRAN_NAME( PSCSTANGRADFACE ,pscstangradface )(
      CHFp_FRA1(gradPhi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(gradBox)
      ,CHFp_INT(gradDir)
      ,CHFp_CONST_REAL(dx) );

#define FORT_PSCSTANGRADFACE FORTRAN_NAME( inlinePSCSTANGRADFACE, inlinePSCSTANGRADFACE)
#define FORTNT_PSCSTANGRADFACE FORTRAN_NAME( PSCSTANGRADFACE, pscstangradface)

inline void FORTRAN_NAME(inlinePSCSTANGRADFACE, inlinePSCSTANGRADFACE)(
      CHFp_FRA1(gradPhi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(gradBox)
      ,CHFp_INT(gradDir)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_PSCSTANGRADFACE");
 FORTRAN_NAME( PSCSTANGRADFACE ,pscstangradface )(
      CHFt_FRA1(gradPhi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_BOX(gradBox)
      ,CHFt_INT(gradDir)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPSCSTANGRADFACE 

#ifndef GUARDFLUXNORMALDIVERGENCE 
#define GUARDFLUXNORMALDIVERGENCE 
// Prototype for Fortran procedure fluxnormaldivergence ...
//
void FORTRAN_NAME( FLUXNORMALDIVERGENCE ,fluxnormaldivergence )(
      CHFp_CONST_FRA(Flux)
      ,CHFp_FRA(div)
      ,CHFp_BOX(gridInt)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_INT(dir) );

#define FORT_FLUXNORMALDIVERGENCE FORTRAN_NAME( inlineFLUXNORMALDIVERGENCE, inlineFLUXNORMALDIVERGENCE)
#define FORTNT_FLUXNORMALDIVERGENCE FORTRAN_NAME( FLUXNORMALDIVERGENCE, fluxnormaldivergence)

inline void FORTRAN_NAME(inlineFLUXNORMALDIVERGENCE, inlineFLUXNORMALDIVERGENCE)(
      CHFp_CONST_FRA(Flux)
      ,CHFp_FRA(div)
      ,CHFp_BOX(gridInt)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_FLUXNORMALDIVERGENCE");
 FORTRAN_NAME( FLUXNORMALDIVERGENCE ,fluxnormaldivergence )(
      CHFt_CONST_FRA(Flux)
      ,CHFt_FRA(div)
      ,CHFt_BOX(gridInt)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_INT(dir) );
}
#endif  // GUARDFLUXNORMALDIVERGENCE 

#ifndef GUARDPSCSDIVERGENCE 
#define GUARDPSCSDIVERGENCE 
// Prototype for Fortran procedure pscsdivergence ...
//
void FORTRAN_NAME( PSCSDIVERGENCE ,pscsdivergence )(
      CHFp_CONST_FRA(uEdge)
      ,CHFp_FRA(div)
      ,CHFp_BOX(gridInt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir) );

#define FORT_PSCSDIVERGENCE FORTRAN_NAME( inlinePSCSDIVERGENCE, inlinePSCSDIVERGENCE)
#define FORTNT_PSCSDIVERGENCE FORTRAN_NAME( PSCSDIVERGENCE, pscsdivergence)

inline void FORTRAN_NAME(inlinePSCSDIVERGENCE, inlinePSCSDIVERGENCE)(
      CHFp_CONST_FRA(uEdge)
      ,CHFp_FRA(div)
      ,CHFp_BOX(gridInt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_PSCSDIVERGENCE");
 FORTRAN_NAME( PSCSDIVERGENCE ,pscsdivergence )(
      CHFt_CONST_FRA(uEdge)
      ,CHFt_FRA(div)
      ,CHFt_BOX(gridInt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(dir) );
}
#endif  // GUARDPSCSDIVERGENCE 

#ifndef GUARDMULT_CFG_NT 
#define GUARDMULT_CFG_NT 
// Prototype for Fortran procedure mult_cfg_nt ...
//
void FORTRAN_NAME( MULT_CFG_NT ,mult_cfg_nt )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(N)
      ,CHFp_FRA(velocity) );

#define FORT_MULT_CFG_NT FORTRAN_NAME( inlineMULT_CFG_NT, inlineMULT_CFG_NT)
#define FORTNT_MULT_CFG_NT FORTRAN_NAME( MULT_CFG_NT, mult_cfg_nt)

inline void FORTRAN_NAME(inlineMULT_CFG_NT, inlineMULT_CFG_NT)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(N)
      ,CHFp_FRA(velocity) )
{
 CH_TIMELEAF("FORT_MULT_CFG_NT");
 FORTRAN_NAME( MULT_CFG_NT ,mult_cfg_nt )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(N)
      ,CHFt_FRA(velocity) );
}
#endif  // GUARDMULT_CFG_NT 

#ifndef GUARDMULT_VEL_NT 
#define GUARDMULT_VEL_NT 
// Prototype for Fortran procedure mult_vel_nt ...
//
void FORTRAN_NAME( MULT_VEL_NT ,mult_vel_nt )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(TwoPirRmaj)
      ,CHFp_FRA(velocity) );

#define FORT_MULT_VEL_NT FORTRAN_NAME( inlineMULT_VEL_NT, inlineMULT_VEL_NT)
#define FORTNT_MULT_VEL_NT FORTRAN_NAME( MULT_VEL_NT, mult_vel_nt)

inline void FORTRAN_NAME(inlineMULT_VEL_NT, inlineMULT_VEL_NT)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(TwoPirRmaj)
      ,CHFp_FRA(velocity) )
{
 CH_TIMELEAF("FORT_MULT_VEL_NT");
 FORTRAN_NAME( MULT_VEL_NT ,mult_vel_nt )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(TwoPirRmaj)
      ,CHFt_FRA(velocity) );
}
#endif  // GUARDMULT_VEL_NT 

#ifndef GUARDMULT_CFG 
#define GUARDMULT_CFG 
// Prototype for Fortran procedure mult_cfg ...
//
void FORTRAN_NAME( MULT_CFG ,mult_cfg )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(factor)
      ,CHFp_FRA(data) );

#define FORT_MULT_CFG FORTRAN_NAME( inlineMULT_CFG, inlineMULT_CFG)
#define FORTNT_MULT_CFG FORTRAN_NAME( MULT_CFG, mult_cfg)

inline void FORTRAN_NAME(inlineMULT_CFG, inlineMULT_CFG)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(factor)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_MULT_CFG");
 FORTRAN_NAME( MULT_CFG ,mult_cfg )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(factor)
      ,CHFt_FRA(data) );
}
#endif  // GUARDMULT_CFG 

#ifndef GUARDMULT_VEL 
#define GUARDMULT_VEL 
// Prototype for Fortran procedure mult_vel ...
//
void FORTRAN_NAME( MULT_VEL ,mult_vel )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(factor)
      ,CHFp_FRA(data) );

#define FORT_MULT_VEL FORTRAN_NAME( inlineMULT_VEL, inlineMULT_VEL)
#define FORTNT_MULT_VEL FORTRAN_NAME( MULT_VEL, mult_vel)

inline void FORTRAN_NAME(inlineMULT_VEL, inlineMULT_VEL)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(factor)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_MULT_VEL");
 FORTRAN_NAME( MULT_VEL ,mult_vel )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(factor)
      ,CHFt_FRA(data) );
}
#endif  // GUARDMULT_VEL 

#ifndef GUARDMULT_RADIAL_COMPONENT 
#define GUARDMULT_RADIAL_COMPONENT 
// Prototype for Fortran procedure mult_radial_component ...
//
void FORTRAN_NAME( MULT_RADIAL_COMPONENT ,mult_radial_component )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(factor)
      ,CHFp_FRA1(data) );

#define FORT_MULT_RADIAL_COMPONENT FORTRAN_NAME( inlineMULT_RADIAL_COMPONENT, inlineMULT_RADIAL_COMPONENT)
#define FORTNT_MULT_RADIAL_COMPONENT FORTRAN_NAME( MULT_RADIAL_COMPONENT, mult_radial_component)

inline void FORTRAN_NAME(inlineMULT_RADIAL_COMPONENT, inlineMULT_RADIAL_COMPONENT)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(factor)
      ,CHFp_FRA1(data) )
{
 CH_TIMELEAF("FORT_MULT_RADIAL_COMPONENT");
 FORTRAN_NAME( MULT_RADIAL_COMPONENT ,mult_radial_component )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(factor)
      ,CHFt_FRA1(data) );
}
#endif  // GUARDMULT_RADIAL_COMPONENT 

#ifndef GUARDCOMPUTE_Q0X 
#define GUARDCOMPUTE_Q0X 
// Prototype for Fortran procedure compute_Q0X ...
//
void FORTRAN_NAME( COMPUTE_Q0X ,compute_q0x )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(integrals)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_REAL(larmor)
      ,CHFp_CONST_INT(include_ExB_drift)
      ,CHFp_CONST_INT(include_magnetic_drift)
      ,CHFp_FRA1(Q) );

#define FORT_COMPUTE_Q0X FORTRAN_NAME( inlineCOMPUTE_Q0X, inlineCOMPUTE_Q0X)
#define FORTNT_COMPUTE_Q0X FORTRAN_NAME( COMPUTE_Q0X, compute_q0x)

inline void FORTRAN_NAME(inlineCOMPUTE_Q0X, inlineCOMPUTE_Q0X)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(integrals)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_REAL(larmor)
      ,CHFp_CONST_INT(include_ExB_drift)
      ,CHFp_CONST_INT(include_magnetic_drift)
      ,CHFp_FRA1(Q) )
{
 CH_TIMELEAF("FORT_COMPUTE_Q0X");
 FORTRAN_NAME( COMPUTE_Q0X ,compute_q0x )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(integrals)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(charge)
      ,CHFt_CONST_REAL(larmor)
      ,CHFt_CONST_INT(include_ExB_drift)
      ,CHFt_CONST_INT(include_magnetic_drift)
      ,CHFt_FRA1(Q) );
}
#endif  // GUARDCOMPUTE_Q0X 

#ifndef GUARDCOMPUTE_Q12 
#define GUARDCOMPUTE_Q12 
// Prototype for Fortran procedure compute_Q12 ...
//
void FORTRAN_NAME( COMPUTE_Q12 ,compute_q12 )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_FRA(nodal_integrals)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_REAL(larmor)
      ,CHFp_CONST_INT(no_par_stream)
      ,CHFp_CONST_INT(include_magnetic_drift)
      ,CHFp_FRA1(Q) );

#define FORT_COMPUTE_Q12 FORTRAN_NAME( inlineCOMPUTE_Q12, inlineCOMPUTE_Q12)
#define FORTNT_COMPUTE_Q12 FORTRAN_NAME( COMPUTE_Q12, compute_q12)

inline void FORTRAN_NAME(inlineCOMPUTE_Q12, inlineCOMPUTE_Q12)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_FRA(nodal_integrals)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_REAL(larmor)
      ,CHFp_CONST_INT(no_par_stream)
      ,CHFp_CONST_INT(include_magnetic_drift)
      ,CHFp_FRA1(Q) )
{
 CH_TIMELEAF("FORT_COMPUTE_Q12");
 FORTRAN_NAME( COMPUTE_Q12 ,compute_q12 )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_FRA(nodal_integrals)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_REAL(charge)
      ,CHFt_CONST_REAL(larmor)
      ,CHFt_CONST_INT(no_par_stream)
      ,CHFt_CONST_INT(include_magnetic_drift)
      ,CHFt_FRA1(Q) );
}
#endif  // GUARDCOMPUTE_Q12 

#ifndef GUARDCOMPUTE_UHAT 
#define GUARDCOMPUTE_UHAT 
// Prototype for Fortran procedure compute_uhat ...
//
void FORTRAN_NAME( COMPUTE_UHAT ,compute_uhat )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(integrals)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_INT(no_zero_order_par)
      ,CHFp_FRA1(uhat) );

#define FORT_COMPUTE_UHAT FORTRAN_NAME( inlineCOMPUTE_UHAT, inlineCOMPUTE_UHAT)
#define FORTNT_COMPUTE_UHAT FORTRAN_NAME( COMPUTE_UHAT, compute_uhat)

inline void FORTRAN_NAME(inlineCOMPUTE_UHAT, inlineCOMPUTE_UHAT)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(integrals)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_INT(no_zero_order_par)
      ,CHFp_FRA1(uhat) )
{
 CH_TIMELEAF("FORT_COMPUTE_UHAT");
 FORTRAN_NAME( COMPUTE_UHAT ,compute_uhat )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(integrals)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_REAL(charge)
      ,CHFt_CONST_INT(no_zero_order_par)
      ,CHFt_FRA1(uhat) );
}
#endif  // GUARDCOMPUTE_UHAT 

#ifndef GUARDCOMPUTE_RADIAL_PROJECTION 
#define GUARDCOMPUTE_RADIAL_PROJECTION 
// Prototype for Fortran procedure compute_radial_projection ...
//
void FORTRAN_NAME( COMPUTE_RADIAL_PROJECTION ,compute_radial_projection )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(vector)
      ,CHFp_CONST_FRA1(vector_r) );

#define FORT_COMPUTE_RADIAL_PROJECTION FORTRAN_NAME( inlineCOMPUTE_RADIAL_PROJECTION, inlineCOMPUTE_RADIAL_PROJECTION)
#define FORTNT_COMPUTE_RADIAL_PROJECTION FORTRAN_NAME( COMPUTE_RADIAL_PROJECTION, compute_radial_projection)

inline void FORTRAN_NAME(inlineCOMPUTE_RADIAL_PROJECTION, inlineCOMPUTE_RADIAL_PROJECTION)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(vector)
      ,CHFp_CONST_FRA1(vector_r) )
{
 CH_TIMELEAF("FORT_COMPUTE_RADIAL_PROJECTION");
 FORTRAN_NAME( COMPUTE_RADIAL_PROJECTION ,compute_radial_projection )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_CONST_FRA(vector)
      ,CHFt_CONST_FRA1(vector_r) );
}
#endif  // GUARDCOMPUTE_RADIAL_PROJECTION 

#ifndef GUARDCOMPUTE_POLOIDAL_PROJECTION 
#define GUARDCOMPUTE_POLOIDAL_PROJECTION 
// Prototype for Fortran procedure compute_poloidal_projection ...
//
void FORTRAN_NAME( COMPUTE_POLOIDAL_PROJECTION ,compute_poloidal_projection )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(vector)
      ,CHFp_CONST_FRA1(vector_pol) );

#define FORT_COMPUTE_POLOIDAL_PROJECTION FORTRAN_NAME( inlineCOMPUTE_POLOIDAL_PROJECTION, inlineCOMPUTE_POLOIDAL_PROJECTION)
#define FORTNT_COMPUTE_POLOIDAL_PROJECTION FORTRAN_NAME( COMPUTE_POLOIDAL_PROJECTION, compute_poloidal_projection)

inline void FORTRAN_NAME(inlineCOMPUTE_POLOIDAL_PROJECTION, inlineCOMPUTE_POLOIDAL_PROJECTION)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(vector)
      ,CHFp_CONST_FRA1(vector_pol) )
{
 CH_TIMELEAF("FORT_COMPUTE_POLOIDAL_PROJECTION");
 FORTRAN_NAME( COMPUTE_POLOIDAL_PROJECTION ,compute_poloidal_projection )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_CONST_FRA(vector)
      ,CHFt_CONST_FRA1(vector_pol) );
}
#endif  // GUARDCOMPUTE_POLOIDAL_PROJECTION 

}

#endif
