#ifndef _DIFFUSIONF_F_H_
#define _DIFFUSIONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_MAPPED_DIFFUSION_COEFFICIENTS 
#define GUARDCOMPUTE_MAPPED_DIFFUSION_COEFFICIENTS 
// Prototype for Fortran procedure compute_mapped_diffusion_coefficients ...
//
void FORTRAN_NAME( COMPUTE_MAPPED_DIFFUSION_COEFFICIENTS ,compute_mapped_diffusion_coefficients )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(n)
      ,CHFp_CONST_FRA(d_tensor)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_FRA(coef) );

#define FORT_COMPUTE_MAPPED_DIFFUSION_COEFFICIENTS FORTRAN_NAME( inlineCOMPUTE_MAPPED_DIFFUSION_COEFFICIENTS, inlineCOMPUTE_MAPPED_DIFFUSION_COEFFICIENTS)
#define FORTNT_COMPUTE_MAPPED_DIFFUSION_COEFFICIENTS FORTRAN_NAME( COMPUTE_MAPPED_DIFFUSION_COEFFICIENTS, compute_mapped_diffusion_coefficients)

inline void FORTRAN_NAME(inlineCOMPUTE_MAPPED_DIFFUSION_COEFFICIENTS, inlineCOMPUTE_MAPPED_DIFFUSION_COEFFICIENTS)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(n)
      ,CHFp_CONST_FRA(d_tensor)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_FRA(coef) )
{
 CH_TIMELEAF("FORT_COMPUTE_MAPPED_DIFFUSION_COEFFICIENTS");
 FORTRAN_NAME( COMPUTE_MAPPED_DIFFUSION_COEFFICIENTS ,compute_mapped_diffusion_coefficients )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(n)
      ,CHFt_CONST_FRA(d_tensor)
      ,CHFt_CONST_FRA(njinverse)
      ,CHFt_FRA(coef) );
}
#endif  // GUARDCOMPUTE_MAPPED_DIFFUSION_COEFFICIENTS 

}

#endif
