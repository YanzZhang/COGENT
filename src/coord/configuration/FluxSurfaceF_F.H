#ifndef _FLUXSURFACEF_F_H_
#define _FLUXSURFACEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDADD_FLUX_SURFACE_ARRAY 
#define GUARDADD_FLUX_SURFACE_ARRAY 
// Prototype for Fortran procedure add_flux_surface_array ...
//
void FORTRAN_NAME( ADD_FLUX_SURFACE_ARRAY ,add_flux_surface_array )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(radial)
      ,CHFp_FRA(data) );

#define FORT_ADD_FLUX_SURFACE_ARRAY FORTRAN_NAME( inlineADD_FLUX_SURFACE_ARRAY, inlineADD_FLUX_SURFACE_ARRAY)
#define FORTNT_ADD_FLUX_SURFACE_ARRAY FORTRAN_NAME( ADD_FLUX_SURFACE_ARRAY, add_flux_surface_array)

inline void FORTRAN_NAME(inlineADD_FLUX_SURFACE_ARRAY, inlineADD_FLUX_SURFACE_ARRAY)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(radial)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_ADD_FLUX_SURFACE_ARRAY");
 FORTRAN_NAME( ADD_FLUX_SURFACE_ARRAY ,add_flux_surface_array )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(radial)
      ,CHFt_FRA(data) );
}
#endif  // GUARDADD_FLUX_SURFACE_ARRAY 

#ifndef GUARDSUBTRACT_FLUX_SURFACE_ARRAY 
#define GUARDSUBTRACT_FLUX_SURFACE_ARRAY 
// Prototype for Fortran procedure subtract_flux_surface_array ...
//
void FORTRAN_NAME( SUBTRACT_FLUX_SURFACE_ARRAY ,subtract_flux_surface_array )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(radial)
      ,CHFp_FRA(data) );

#define FORT_SUBTRACT_FLUX_SURFACE_ARRAY FORTRAN_NAME( inlineSUBTRACT_FLUX_SURFACE_ARRAY, inlineSUBTRACT_FLUX_SURFACE_ARRAY)
#define FORTNT_SUBTRACT_FLUX_SURFACE_ARRAY FORTRAN_NAME( SUBTRACT_FLUX_SURFACE_ARRAY, subtract_flux_surface_array)

inline void FORTRAN_NAME(inlineSUBTRACT_FLUX_SURFACE_ARRAY, inlineSUBTRACT_FLUX_SURFACE_ARRAY)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(radial)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_SUBTRACT_FLUX_SURFACE_ARRAY");
 FORTRAN_NAME( SUBTRACT_FLUX_SURFACE_ARRAY ,subtract_flux_surface_array )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(radial)
      ,CHFt_FRA(data) );
}
#endif  // GUARDSUBTRACT_FLUX_SURFACE_ARRAY 

#ifndef GUARDMULTIPLY_FLUX_SURFACE_ARRAY 
#define GUARDMULTIPLY_FLUX_SURFACE_ARRAY 
// Prototype for Fortran procedure multiply_flux_surface_array ...
//
void FORTRAN_NAME( MULTIPLY_FLUX_SURFACE_ARRAY ,multiply_flux_surface_array )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(radial)
      ,CHFp_FRA(data) );

#define FORT_MULTIPLY_FLUX_SURFACE_ARRAY FORTRAN_NAME( inlineMULTIPLY_FLUX_SURFACE_ARRAY, inlineMULTIPLY_FLUX_SURFACE_ARRAY)
#define FORTNT_MULTIPLY_FLUX_SURFACE_ARRAY FORTRAN_NAME( MULTIPLY_FLUX_SURFACE_ARRAY, multiply_flux_surface_array)

inline void FORTRAN_NAME(inlineMULTIPLY_FLUX_SURFACE_ARRAY, inlineMULTIPLY_FLUX_SURFACE_ARRAY)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(radial)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_MULTIPLY_FLUX_SURFACE_ARRAY");
 FORTRAN_NAME( MULTIPLY_FLUX_SURFACE_ARRAY ,multiply_flux_surface_array )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(radial)
      ,CHFt_FRA(data) );
}
#endif  // GUARDMULTIPLY_FLUX_SURFACE_ARRAY 

#ifndef GUARDDIVIDE_FLUX_SURFACE_ARRAY 
#define GUARDDIVIDE_FLUX_SURFACE_ARRAY 
// Prototype for Fortran procedure divide_flux_surface_array ...
//
void FORTRAN_NAME( DIVIDE_FLUX_SURFACE_ARRAY ,divide_flux_surface_array )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(radial)
      ,CHFp_FRA(data) );

#define FORT_DIVIDE_FLUX_SURFACE_ARRAY FORTRAN_NAME( inlineDIVIDE_FLUX_SURFACE_ARRAY, inlineDIVIDE_FLUX_SURFACE_ARRAY)
#define FORTNT_DIVIDE_FLUX_SURFACE_ARRAY FORTRAN_NAME( DIVIDE_FLUX_SURFACE_ARRAY, divide_flux_surface_array)

inline void FORTRAN_NAME(inlineDIVIDE_FLUX_SURFACE_ARRAY, inlineDIVIDE_FLUX_SURFACE_ARRAY)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(radial)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_DIVIDE_FLUX_SURFACE_ARRAY");
 FORTRAN_NAME( DIVIDE_FLUX_SURFACE_ARRAY ,divide_flux_surface_array )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(radial)
      ,CHFt_FRA(data) );
}
#endif  // GUARDDIVIDE_FLUX_SURFACE_ARRAY 

}

#endif
