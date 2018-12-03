#ifndef _MBSOLVERUTILF_F_H_
#define _MBSOLVERUTILF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_MAPPED_COEFFICIENTS 
#define GUARDCOMPUTE_MAPPED_COEFFICIENTS 
// Prototype for Fortran procedure compute_mapped_coefficients ...
//
void FORTRAN_NAME( COMPUTE_MAPPED_COEFFICIENTS ,compute_mapped_coefficients )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(unmapped_coef)
      ,CHFp_CONST_FRA(n)
      ,CHFp_CONST_FRA1(jinverse)
      ,CHFp_FRA(coef) );

#define FORT_COMPUTE_MAPPED_COEFFICIENTS FORTRAN_NAME( inlineCOMPUTE_MAPPED_COEFFICIENTS, inlineCOMPUTE_MAPPED_COEFFICIENTS)
#define FORTNT_COMPUTE_MAPPED_COEFFICIENTS FORTRAN_NAME( COMPUTE_MAPPED_COEFFICIENTS, compute_mapped_coefficients)

inline void FORTRAN_NAME(inlineCOMPUTE_MAPPED_COEFFICIENTS, inlineCOMPUTE_MAPPED_COEFFICIENTS)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(unmapped_coef)
      ,CHFp_CONST_FRA(n)
      ,CHFp_CONST_FRA1(jinverse)
      ,CHFp_FRA(coef) )
{
 CH_TIMELEAF("FORT_COMPUTE_MAPPED_COEFFICIENTS");
 FORTRAN_NAME( COMPUTE_MAPPED_COEFFICIENTS ,compute_mapped_coefficients )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(unmapped_coef)
      ,CHFt_CONST_FRA(n)
      ,CHFt_CONST_FRA1(jinverse)
      ,CHFt_FRA(coef) );
}
#endif  // GUARDCOMPUTE_MAPPED_COEFFICIENTS 

}

#endif
