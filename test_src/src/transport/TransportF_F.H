#ifndef _TRANSPORTF_F_H_
#define _TRANSPORTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEVAL_ANOM_FLUX 
#define GUARDEVAL_ANOM_FLUX 
// Prototype for Fortran procedure eval_anom_flux ...
//
void FORTRAN_NAME( EVAL_ANOM_FLUX ,eval_anom_flux )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(lame_face)
      ,CHFp_CONST_FRA(D_cfg)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_INT(Nr)
      ,CHFp_CONST_INT(simpleDiffusion)
      ,CHFp_CONST_FRA1(fB)
      ,CHFp_CONST_FRA1(dfBdmu)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_FRA1(N)
      ,CHFp_CONST_FRA1(U)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(C)
      ,CHFp_CONST_FRA1(P)
      ,CHFp_FRA1(flux) );

#define FORT_EVAL_ANOM_FLUX FORTRAN_NAME( inlineEVAL_ANOM_FLUX, inlineEVAL_ANOM_FLUX)
#define FORTNT_EVAL_ANOM_FLUX FORTRAN_NAME( EVAL_ANOM_FLUX, eval_anom_flux)

inline void FORTRAN_NAME(inlineEVAL_ANOM_FLUX, inlineEVAL_ANOM_FLUX)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(lame_face)
      ,CHFp_CONST_FRA(D_cfg)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_INT(Nr)
      ,CHFp_CONST_INT(simpleDiffusion)
      ,CHFp_CONST_FRA1(fB)
      ,CHFp_CONST_FRA1(dfBdmu)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_FRA1(N)
      ,CHFp_CONST_FRA1(U)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(C)
      ,CHFp_CONST_FRA1(P)
      ,CHFp_FRA1(flux) )
{
 CH_TIMELEAF("FORT_EVAL_ANOM_FLUX");
 FORTRAN_NAME( EVAL_ANOM_FLUX ,eval_anom_flux )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_FRA(lame_face)
      ,CHFt_CONST_FRA(D_cfg)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_INT(Nr)
      ,CHFt_CONST_INT(simpleDiffusion)
      ,CHFt_CONST_FRA1(fB)
      ,CHFt_CONST_FRA1(dfBdmu)
      ,CHFt_CONST_FRA1(B)
      ,CHFt_CONST_FRA1(N)
      ,CHFt_CONST_FRA1(U)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_CONST_FRA1(C)
      ,CHFt_CONST_FRA1(P)
      ,CHFt_FRA1(flux) );
}
#endif  // GUARDEVAL_ANOM_FLUX 

#ifndef GUARDEVAL_ANOM_FLUX_NOT_ALIGNED 
#define GUARDEVAL_ANOM_FLUX_NOT_ALIGNED 
// Prototype for Fortran procedure eval_anom_flux_not_aligned ...
//
void FORTRAN_NAME( EVAL_ANOM_FLUX_NOT_ALIGNED ,eval_anom_flux_not_aligned )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(NJinv)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(D_cfg)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_INT(Nr)
      ,CHFp_CONST_INT(simpleDiffusion)
      ,CHFp_CONST_FRA1(fB)
      ,CHFp_CONST_FRA1(dfBdmu)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_FRA1(N)
      ,CHFp_CONST_FRA1(U)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(C)
      ,CHFp_CONST_FRA1(P)
      ,CHFp_FRA(flux) );

#define FORT_EVAL_ANOM_FLUX_NOT_ALIGNED FORTRAN_NAME( inlineEVAL_ANOM_FLUX_NOT_ALIGNED, inlineEVAL_ANOM_FLUX_NOT_ALIGNED)
#define FORTNT_EVAL_ANOM_FLUX_NOT_ALIGNED FORTRAN_NAME( EVAL_ANOM_FLUX_NOT_ALIGNED, eval_anom_flux_not_aligned)

inline void FORTRAN_NAME(inlineEVAL_ANOM_FLUX_NOT_ALIGNED, inlineEVAL_ANOM_FLUX_NOT_ALIGNED)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(NJinv)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(D_cfg)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_INT(Nr)
      ,CHFp_CONST_INT(simpleDiffusion)
      ,CHFp_CONST_FRA1(fB)
      ,CHFp_CONST_FRA1(dfBdmu)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_FRA1(N)
      ,CHFp_CONST_FRA1(U)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(C)
      ,CHFp_CONST_FRA1(P)
      ,CHFp_FRA(flux) )
{
 CH_TIMELEAF("FORT_EVAL_ANOM_FLUX_NOT_ALIGNED");
 FORTRAN_NAME( EVAL_ANOM_FLUX_NOT_ALIGNED ,eval_anom_flux_not_aligned )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_FRA(NJinv)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_CONST_FRA(D_cfg)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_INT(Nr)
      ,CHFt_CONST_INT(simpleDiffusion)
      ,CHFt_CONST_FRA1(fB)
      ,CHFt_CONST_FRA1(dfBdmu)
      ,CHFt_CONST_FRA1(B)
      ,CHFt_CONST_FRA1(N)
      ,CHFt_CONST_FRA1(U)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_CONST_FRA1(C)
      ,CHFt_CONST_FRA1(P)
      ,CHFt_FRA(flux) );
}
#endif  // GUARDEVAL_ANOM_FLUX_NOT_ALIGNED 

#ifndef GUARDDFBDMU_CELL_CENTER 
#define GUARDDFBDMU_CELL_CENTER 
// Prototype for Fortran procedure dfBdmu_cell_center ...
//
void FORTRAN_NAME( DFBDMU_CELL_CENTER ,dfbdmu_cell_center )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(nmu)
      ,CHFp_CONST_FRA1(fB)
      ,CHFp_FRA1(dfB_dmu) );

#define FORT_DFBDMU_CELL_CENTER FORTRAN_NAME( inlineDFBDMU_CELL_CENTER, inlineDFBDMU_CELL_CENTER)
#define FORTNT_DFBDMU_CELL_CENTER FORTRAN_NAME( DFBDMU_CELL_CENTER, dfbdmu_cell_center)

inline void FORTRAN_NAME(inlineDFBDMU_CELL_CENTER, inlineDFBDMU_CELL_CENTER)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(nmu)
      ,CHFp_CONST_FRA1(fB)
      ,CHFp_FRA1(dfB_dmu) )
{
 CH_TIMELEAF("FORT_DFBDMU_CELL_CENTER");
 FORTRAN_NAME( DFBDMU_CELL_CENTER ,dfbdmu_cell_center )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(nmu)
      ,CHFt_CONST_FRA1(fB)
      ,CHFt_FRA1(dfB_dmu) );
}
#endif  // GUARDDFBDMU_CELL_CENTER 

#ifndef GUARDLAME_COEFFICIENTS 
#define GUARDLAME_COEFFICIENTS 
// Prototype for Fortran procedure lame_coefficients ...
//
void FORTRAN_NAME( LAME_COEFFICIENTS ,lame_coefficients )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(N)
      ,CHFp_CONST_FRA(dXdq)
      ,CHFp_FRA(Lame_coeff) );

#define FORT_LAME_COEFFICIENTS FORTRAN_NAME( inlineLAME_COEFFICIENTS, inlineLAME_COEFFICIENTS)
#define FORTNT_LAME_COEFFICIENTS FORTRAN_NAME( LAME_COEFFICIENTS, lame_coefficients)

inline void FORTRAN_NAME(inlineLAME_COEFFICIENTS, inlineLAME_COEFFICIENTS)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(N)
      ,CHFp_CONST_FRA(dXdq)
      ,CHFp_FRA(Lame_coeff) )
{
 CH_TIMELEAF("FORT_LAME_COEFFICIENTS");
 FORTRAN_NAME( LAME_COEFFICIENTS ,lame_coefficients )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(N)
      ,CHFt_CONST_FRA(dXdq)
      ,CHFt_FRA(Lame_coeff) );
}
#endif  // GUARDLAME_COEFFICIENTS 

#ifndef GUARDEVAL_BETA 
#define GUARDEVAL_BETA 
// Prototype for Fortran procedure eval_beta ...
//
void FORTRAN_NAME( EVAL_BETA ,eval_beta )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_VR(D_kinet)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_INT(Nr)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_FRA1(N)
      ,CHFp_CONST_FRA1(U)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(C)
      ,CHFp_CONST_FRA1(P)
      ,CHFp_FRA1(beta) );

#define FORT_EVAL_BETA FORTRAN_NAME( inlineEVAL_BETA, inlineEVAL_BETA)
#define FORTNT_EVAL_BETA FORTRAN_NAME( EVAL_BETA, eval_beta)

inline void FORTRAN_NAME(inlineEVAL_BETA, inlineEVAL_BETA)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_VR(D_kinet)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_INT(Nr)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_FRA1(N)
      ,CHFp_CONST_FRA1(U)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(C)
      ,CHFp_CONST_FRA1(P)
      ,CHFp_FRA1(beta) )
{
 CH_TIMELEAF("FORT_EVAL_BETA");
 FORTRAN_NAME( EVAL_BETA ,eval_beta )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_VR(D_kinet)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_INT(Nr)
      ,CHFt_CONST_FRA1(B)
      ,CHFt_CONST_FRA1(N)
      ,CHFt_CONST_FRA1(U)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_CONST_FRA1(C)
      ,CHFt_CONST_FRA1(P)
      ,CHFt_FRA1(beta) );
}
#endif  // GUARDEVAL_BETA 

}

#endif
