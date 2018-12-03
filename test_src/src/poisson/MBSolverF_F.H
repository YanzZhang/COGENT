#ifndef _MBSOLVERF_F_H_
#define _MBSOLVERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDACCUM_FLUX_STENCIL4 
#define GUARDACCUM_FLUX_STENCIL4 
// Prototype for Fortran procedure accum_flux_stencil4 ...
//
void FORTRAN_NAME( ACCUM_FLUX_STENCIL4 ,accum_flux_stencil4 )(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(deriv_dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_FRA1(coef)
      ,CHFp_CONST_INTVECT(global)
      ,CHFp_FRA1(sum) );

#define FORT_ACCUM_FLUX_STENCIL4 FORTRAN_NAME( inlineACCUM_FLUX_STENCIL4, inlineACCUM_FLUX_STENCIL4)
#define FORTNT_ACCUM_FLUX_STENCIL4 FORTRAN_NAME( ACCUM_FLUX_STENCIL4, accum_flux_stencil4)

inline void FORTRAN_NAME(inlineACCUM_FLUX_STENCIL4, inlineACCUM_FLUX_STENCIL4)(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(deriv_dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_FRA1(coef)
      ,CHFp_CONST_INTVECT(global)
      ,CHFp_FRA1(sum) )
{
 CH_TIMELEAF("FORT_ACCUM_FLUX_STENCIL4");
 FORTRAN_NAME( ACCUM_FLUX_STENCIL4 ,accum_flux_stencil4 )(
      CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(deriv_dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_FRA1(coef)
      ,CHFt_CONST_INTVECT(global)
      ,CHFt_FRA1(sum) );
}
#endif  // GUARDACCUM_FLUX_STENCIL4 

#ifndef GUARDACCUM_FLUX_STENCIL2 
#define GUARDACCUM_FLUX_STENCIL2 
// Prototype for Fortran procedure accum_flux_stencil2 ...
//
void FORTRAN_NAME( ACCUM_FLUX_STENCIL2 ,accum_flux_stencil2 )(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(deriv_dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_FRA1(coef)
      ,CHFp_CONST_INTVECT(global)
      ,CHFp_FRA1(sum) );

#define FORT_ACCUM_FLUX_STENCIL2 FORTRAN_NAME( inlineACCUM_FLUX_STENCIL2, inlineACCUM_FLUX_STENCIL2)
#define FORTNT_ACCUM_FLUX_STENCIL2 FORTRAN_NAME( ACCUM_FLUX_STENCIL2, accum_flux_stencil2)

inline void FORTRAN_NAME(inlineACCUM_FLUX_STENCIL2, inlineACCUM_FLUX_STENCIL2)(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(deriv_dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_FRA1(coef)
      ,CHFp_CONST_INTVECT(global)
      ,CHFp_FRA1(sum) )
{
 CH_TIMELEAF("FORT_ACCUM_FLUX_STENCIL2");
 FORTRAN_NAME( ACCUM_FLUX_STENCIL2 ,accum_flux_stencil2 )(
      CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(deriv_dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_FRA1(coef)
      ,CHFt_CONST_INTVECT(global)
      ,CHFt_FRA1(sum) );
}
#endif  // GUARDACCUM_FLUX_STENCIL2 

}

#endif
