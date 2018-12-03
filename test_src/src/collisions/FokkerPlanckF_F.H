#ifndef _FOKKERPLANCKF_F_H_
#define _FOKKERPLANCKF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEVAL_DV 
#define GUARDEVAL_DV 
// Prototype for Fortran procedure eval_dv ...
//
void FORTRAN_NAME( EVAL_DV ,eval_dv )(
      CHFp_FRA1(Dv)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) );

#define FORT_EVAL_DV FORTRAN_NAME( inlineEVAL_DV, inlineEVAL_DV)
#define FORTNT_EVAL_DV FORTRAN_NAME( EVAL_DV, eval_dv)

inline void FORTRAN_NAME(inlineEVAL_DV, inlineEVAL_DV)(
      CHFp_FRA1(Dv)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) )
{
 CH_TIMELEAF("FORT_EVAL_DV");
 FORTRAN_NAME( EVAL_DV ,eval_dv )(
      CHFt_FRA1(Dv)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(mTP)
      ,CHFt_CONST_REAL(mFP)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu) );
}
#endif  // GUARDEVAL_DV 

#ifndef GUARDEVAL_DMU 
#define GUARDEVAL_DMU 
// Prototype for Fortran procedure eval_dmu ...
//
void FORTRAN_NAME( EVAL_DMU ,eval_dmu )(
      CHFp_FRA1(Dmu)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) );

#define FORT_EVAL_DMU FORTRAN_NAME( inlineEVAL_DMU, inlineEVAL_DMU)
#define FORTNT_EVAL_DMU FORTRAN_NAME( EVAL_DMU, eval_dmu)

inline void FORTRAN_NAME(inlineEVAL_DMU, inlineEVAL_DMU)(
      CHFp_FRA1(Dmu)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) )
{
 CH_TIMELEAF("FORT_EVAL_DMU");
 FORTRAN_NAME( EVAL_DMU ,eval_dmu )(
      CHFt_FRA1(Dmu)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(mTP)
      ,CHFt_CONST_REAL(mFP)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu) );
}
#endif  // GUARDEVAL_DMU 

#ifndef GUARDEVAL_DVV 
#define GUARDEVAL_DVV 
// Prototype for Fortran procedure eval_dvv ...
//
void FORTRAN_NAME( EVAL_DVV ,eval_dvv )(
      CHFp_FRA1(Dvv)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) );

#define FORT_EVAL_DVV FORTRAN_NAME( inlineEVAL_DVV, inlineEVAL_DVV)
#define FORTNT_EVAL_DVV FORTRAN_NAME( EVAL_DVV, eval_dvv)

inline void FORTRAN_NAME(inlineEVAL_DVV, inlineEVAL_DVV)(
      CHFp_FRA1(Dvv)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) )
{
 CH_TIMELEAF("FORT_EVAL_DVV");
 FORTRAN_NAME( EVAL_DVV ,eval_dvv )(
      CHFt_FRA1(Dvv)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(mTP)
      ,CHFt_CONST_REAL(mFP)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu) );
}
#endif  // GUARDEVAL_DVV 

#ifndef GUARDEVAL_DMUMU 
#define GUARDEVAL_DMUMU 
// Prototype for Fortran procedure eval_dmumu ...
//
void FORTRAN_NAME( EVAL_DMUMU ,eval_dmumu )(
      CHFp_FRA1(Dmumu)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) );

#define FORT_EVAL_DMUMU FORTRAN_NAME( inlineEVAL_DMUMU, inlineEVAL_DMUMU)
#define FORTNT_EVAL_DMUMU FORTRAN_NAME( EVAL_DMUMU, eval_dmumu)

inline void FORTRAN_NAME(inlineEVAL_DMUMU, inlineEVAL_DMUMU)(
      CHFp_FRA1(Dmumu)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) )
{
 CH_TIMELEAF("FORT_EVAL_DMUMU");
 FORTRAN_NAME( EVAL_DMUMU ,eval_dmumu )(
      CHFt_FRA1(Dmumu)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(mTP)
      ,CHFt_CONST_REAL(mFP)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu) );
}
#endif  // GUARDEVAL_DMUMU 

#ifndef GUARDEVAL_DVMU 
#define GUARDEVAL_DVMU 
// Prototype for Fortran procedure eval_dvmu ...
//
void FORTRAN_NAME( EVAL_DVMU ,eval_dvmu )(
      CHFp_FRA1(Dvmu)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) );

#define FORT_EVAL_DVMU FORTRAN_NAME( inlineEVAL_DVMU, inlineEVAL_DVMU)
#define FORTNT_EVAL_DVMU FORTRAN_NAME( EVAL_DVMU, eval_dvmu)

inline void FORTRAN_NAME(inlineEVAL_DVMU, inlineEVAL_DVMU)(
      CHFp_FRA1(Dvmu)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(mTP)
      ,CHFp_CONST_REAL(mFP)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) )
{
 CH_TIMELEAF("FORT_EVAL_DVMU");
 FORTRAN_NAME( EVAL_DVMU ,eval_dvmu )(
      CHFt_FRA1(Dvmu)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(mTP)
      ,CHFt_CONST_REAL(mFP)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu) );
}
#endif  // GUARDEVAL_DVMU 

#ifndef GUARDFLUX_VPAR 
#define GUARDFLUX_VPAR 
// Prototype for Fortran procedure flux_vpar ...
//
void FORTRAN_NAME( FLUX_VPAR ,flux_vpar )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(Dv)
      ,CHFp_CONST_FRA1(Dvv)
      ,CHFp_CONST_FRA1(Dvmu)
      ,CHFp_CONST_FRA(dfn)
      ,CHFp_CONST_FRA(dfn_f)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_INT(flag) );

#define FORT_FLUX_VPAR FORTRAN_NAME( inlineFLUX_VPAR, inlineFLUX_VPAR)
#define FORTNT_FLUX_VPAR FORTRAN_NAME( FLUX_VPAR, flux_vpar)

inline void FORTRAN_NAME(inlineFLUX_VPAR, inlineFLUX_VPAR)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(Dv)
      ,CHFp_CONST_FRA1(Dvv)
      ,CHFp_CONST_FRA1(Dvmu)
      ,CHFp_CONST_FRA(dfn)
      ,CHFp_CONST_FRA(dfn_f)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_INT(flag) )
{
 CH_TIMELEAF("FORT_FLUX_VPAR");
 FORTRAN_NAME( FLUX_VPAR ,flux_vpar )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA1(Dv)
      ,CHFt_CONST_FRA1(Dvv)
      ,CHFt_CONST_FRA1(Dvmu)
      ,CHFt_CONST_FRA(dfn)
      ,CHFt_CONST_FRA(dfn_f)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu)
      ,CHFt_CONST_INT(flag) );
}
#endif  // GUARDFLUX_VPAR 

#ifndef GUARDFLUX_MU 
#define GUARDFLUX_MU 
// Prototype for Fortran procedure flux_mu ...
//
void FORTRAN_NAME( FLUX_MU ,flux_mu )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(Dmu)
      ,CHFp_CONST_FRA1(Dmumu)
      ,CHFp_CONST_FRA1(Dvmu)
      ,CHFp_CONST_FRA(dfn)
      ,CHFp_CONST_FRA(dfn_f)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_INT(flag) );

#define FORT_FLUX_MU FORTRAN_NAME( inlineFLUX_MU, inlineFLUX_MU)
#define FORTNT_FLUX_MU FORTRAN_NAME( FLUX_MU, flux_mu)

inline void FORTRAN_NAME(inlineFLUX_MU, inlineFLUX_MU)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(Dmu)
      ,CHFp_CONST_FRA1(Dmumu)
      ,CHFp_CONST_FRA1(Dvmu)
      ,CHFp_CONST_FRA(dfn)
      ,CHFp_CONST_FRA(dfn_f)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_INT(flag) )
{
 CH_TIMELEAF("FORT_FLUX_MU");
 FORTRAN_NAME( FLUX_MU ,flux_mu )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA1(Dmu)
      ,CHFt_CONST_FRA1(Dmumu)
      ,CHFt_CONST_FRA1(Dvmu)
      ,CHFt_CONST_FRA(dfn)
      ,CHFt_CONST_FRA(dfn_f)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu)
      ,CHFt_CONST_INT(flag) );
}
#endif  // GUARDFLUX_MU 

#ifndef GUARDEVAL_RHS_COMPONENTS 
#define GUARDEVAL_RHS_COMPONENTS 
// Prototype for Fortran procedure eval_rhs_components ...
//
void FORTRAN_NAME( EVAL_RHS_COMPONENTS ,eval_rhs_components )(
      CHFp_FRA(rhs)
      ,CHFp_CONST_FRA(flux_vpar)
      ,CHFp_CONST_FRA(flux_mu)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx) );

#define FORT_EVAL_RHS_COMPONENTS FORTRAN_NAME( inlineEVAL_RHS_COMPONENTS, inlineEVAL_RHS_COMPONENTS)
#define FORTNT_EVAL_RHS_COMPONENTS FORTRAN_NAME( EVAL_RHS_COMPONENTS, eval_rhs_components)

inline void FORTRAN_NAME(inlineEVAL_RHS_COMPONENTS, inlineEVAL_RHS_COMPONENTS)(
      CHFp_FRA(rhs)
      ,CHFp_CONST_FRA(flux_vpar)
      ,CHFp_CONST_FRA(flux_mu)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx) )
{
 CH_TIMELEAF("FORT_EVAL_RHS_COMPONENTS");
 FORTRAN_NAME( EVAL_RHS_COMPONENTS ,eval_rhs_components )(
      CHFt_FRA(rhs)
      ,CHFt_CONST_FRA(flux_vpar)
      ,CHFt_CONST_FRA(flux_mu)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx) );
}
#endif  // GUARDEVAL_RHS_COMPONENTS 

#ifndef GUARDEVAL_RHS 
#define GUARDEVAL_RHS 
// Prototype for Fortran procedure eval_rhs ...
//
void FORTRAN_NAME( EVAL_RHS ,eval_rhs )(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA(flux_vpar)
      ,CHFp_CONST_FRA(flux_mu)
      ,CHFp_CONST_FRA1(e)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(nu) );

#define FORT_EVAL_RHS FORTRAN_NAME( inlineEVAL_RHS, inlineEVAL_RHS)
#define FORTNT_EVAL_RHS FORTRAN_NAME( EVAL_RHS, eval_rhs)

inline void FORTRAN_NAME(inlineEVAL_RHS, inlineEVAL_RHS)(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA(flux_vpar)
      ,CHFp_CONST_FRA(flux_mu)
      ,CHFp_CONST_FRA1(e)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(nu) )
{
 CH_TIMELEAF("FORT_EVAL_RHS");
 FORTRAN_NAME( EVAL_RHS ,eval_rhs )(
      CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA(flux_vpar)
      ,CHFt_CONST_FRA(flux_mu)
      ,CHFt_CONST_FRA1(e)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu)
      ,CHFt_CONST_REAL(nu) );
}
#endif  // GUARDEVAL_RHS 

}

#endif
