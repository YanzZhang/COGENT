#ifndef _GKPOISSONF_F_H_
#define _GKPOISSONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_CONDUCTIVITY 
#define GUARDCOMPUTE_CONDUCTIVITY 
// Prototype for Fortran procedure compute_conductivity ...
//
void FORTRAN_NAME( COMPUTE_CONDUCTIVITY ,compute_conductivity )(
      CHFp_BOX(box)
      ,CHFp_FRA1(sigma)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_REAL(sigma_max) );

#define FORT_COMPUTE_CONDUCTIVITY FORTRAN_NAME( inlineCOMPUTE_CONDUCTIVITY, inlineCOMPUTE_CONDUCTIVITY)
#define FORTNT_COMPUTE_CONDUCTIVITY FORTRAN_NAME( COMPUTE_CONDUCTIVITY, compute_conductivity)

inline void FORTRAN_NAME(inlineCOMPUTE_CONDUCTIVITY, inlineCOMPUTE_CONDUCTIVITY)(
      CHFp_BOX(box)
      ,CHFp_FRA1(sigma)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_REAL(sigma_max) )
{
 CH_TIMELEAF("FORT_COMPUTE_CONDUCTIVITY");
 FORTRAN_NAME( COMPUTE_CONDUCTIVITY ,compute_conductivity )(
      CHFt_BOX(box)
      ,CHFt_FRA1(sigma)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_CONST_REAL(coeff)
      ,CHFt_CONST_REAL(sigma_max) );
}
#endif  // GUARDCOMPUTE_CONDUCTIVITY 

}

#endif
