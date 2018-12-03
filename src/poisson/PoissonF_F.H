#ifndef _POISSONF_F_H_
#define _POISSONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_LAPLACIAN_COEFFICIENTS 
#define GUARDCOMPUTE_LAPLACIAN_COEFFICIENTS 
// Prototype for Fortran procedure compute_laplacian_coefficients ...
//
void FORTRAN_NAME( COMPUTE_LAPLACIAN_COEFFICIENTS ,compute_laplacian_coefficients )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(n)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_FRA(coef) );

#define FORT_COMPUTE_LAPLACIAN_COEFFICIENTS FORTRAN_NAME( inlineCOMPUTE_LAPLACIAN_COEFFICIENTS, inlineCOMPUTE_LAPLACIAN_COEFFICIENTS)
#define FORTNT_COMPUTE_LAPLACIAN_COEFFICIENTS FORTRAN_NAME( COMPUTE_LAPLACIAN_COEFFICIENTS, compute_laplacian_coefficients)

inline void FORTRAN_NAME(inlineCOMPUTE_LAPLACIAN_COEFFICIENTS, inlineCOMPUTE_LAPLACIAN_COEFFICIENTS)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(n)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_FRA(coef) )
{
 CH_TIMELEAF("FORT_COMPUTE_LAPLACIAN_COEFFICIENTS");
 FORTRAN_NAME( COMPUTE_LAPLACIAN_COEFFICIENTS ,compute_laplacian_coefficients )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(n)
      ,CHFt_CONST_FRA(njinverse)
      ,CHFt_FRA(coef) );
}
#endif  // GUARDCOMPUTE_LAPLACIAN_COEFFICIENTS 

}

#endif
