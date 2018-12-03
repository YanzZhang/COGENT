#ifndef _COLLISIONSF_F_H_
#define _COLLISIONSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEVALUATE_TP_LORENTZ_CONST_NUD 
#define GUARDEVALUATE_TP_LORENTZ_CONST_NUD 
// Prototype for Fortran procedure evaluate_tp_lorentz_const_nud ...
//
void FORTRAN_NAME( EVALUATE_TP_LORENTZ_CONST_NUD ,evaluate_tp_lorentz_const_nud )(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m) );

#define FORT_EVALUATE_TP_LORENTZ_CONST_NUD FORTRAN_NAME( inlineEVALUATE_TP_LORENTZ_CONST_NUD, inlineEVALUATE_TP_LORENTZ_CONST_NUD)
#define FORTNT_EVALUATE_TP_LORENTZ_CONST_NUD FORTRAN_NAME( EVALUATE_TP_LORENTZ_CONST_NUD, evaluate_tp_lorentz_const_nud)

inline void FORTRAN_NAME(inlineEVALUATE_TP_LORENTZ_CONST_NUD, inlineEVALUATE_TP_LORENTZ_CONST_NUD)(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m) )
{
 CH_TIMELEAF("FORT_EVALUATE_TP_LORENTZ_CONST_NUD");
 FORTRAN_NAME( EVALUATE_TP_LORENTZ_CONST_NUD ,evaluate_tp_lorentz_const_nud )(
      CHFt_FRA(ClsFlux)
      ,CHFt_CONST_FRA1(f)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(m) );
}
#endif  // GUARDEVALUATE_TP_LORENTZ_CONST_NUD 

#ifndef GUARDEVALUATE_TP_LORENTZ 
#define GUARDEVALUATE_TP_LORENTZ 
// Prototype for Fortran procedure evaluate_tp_lorentz ...
//
void FORTRAN_NAME( EVALUATE_TP_LORENTZ ,evaluate_tp_lorentz )(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m) );

#define FORT_EVALUATE_TP_LORENTZ FORTRAN_NAME( inlineEVALUATE_TP_LORENTZ, inlineEVALUATE_TP_LORENTZ)
#define FORTNT_EVALUATE_TP_LORENTZ FORTRAN_NAME( EVALUATE_TP_LORENTZ, evaluate_tp_lorentz)

inline void FORTRAN_NAME(inlineEVALUATE_TP_LORENTZ, inlineEVALUATE_TP_LORENTZ)(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m) )
{
 CH_TIMELEAF("FORT_EVALUATE_TP_LORENTZ");
 FORTRAN_NAME( EVALUATE_TP_LORENTZ ,evaluate_tp_lorentz )(
      CHFt_FRA(ClsFlux)
      ,CHFt_CONST_FRA1(f)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(m) );
}
#endif  // GUARDEVALUATE_TP_LORENTZ 

#ifndef GUARDEVALUATE_TP_LORENTZ_SPECIES 
#define GUARDEVALUATE_TP_LORENTZ_SPECIES 
// Prototype for Fortran procedure evaluate_tp_lorentz_species ...
//
void FORTRAN_NAME( EVALUATE_TP_LORENTZ_SPECIES ,evaluate_tp_lorentz_species )(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m)
      ,CHFp_CONST_REAL(mb) );

#define FORT_EVALUATE_TP_LORENTZ_SPECIES FORTRAN_NAME( inlineEVALUATE_TP_LORENTZ_SPECIES, inlineEVALUATE_TP_LORENTZ_SPECIES)
#define FORTNT_EVALUATE_TP_LORENTZ_SPECIES FORTRAN_NAME( EVALUATE_TP_LORENTZ_SPECIES, evaluate_tp_lorentz_species)

inline void FORTRAN_NAME(inlineEVALUATE_TP_LORENTZ_SPECIES, inlineEVALUATE_TP_LORENTZ_SPECIES)(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m)
      ,CHFp_CONST_REAL(mb) )
{
 CH_TIMELEAF("FORT_EVALUATE_TP_LORENTZ_SPECIES");
 FORTRAN_NAME( EVALUATE_TP_LORENTZ_SPECIES ,evaluate_tp_lorentz_species )(
      CHFt_FRA(ClsFlux)
      ,CHFt_CONST_FRA1(f)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(m)
      ,CHFt_CONST_REAL(mb) );
}
#endif  // GUARDEVALUATE_TP_LORENTZ_SPECIES 

#ifndef GUARDEVALUATE_TP_ENERG_DIFF 
#define GUARDEVALUATE_TP_ENERG_DIFF 
// Prototype for Fortran procedure evaluate_tp_energ_diff ...
//
void FORTRAN_NAME( EVALUATE_TP_ENERG_DIFF ,evaluate_tp_energ_diff )(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m) );

#define FORT_EVALUATE_TP_ENERG_DIFF FORTRAN_NAME( inlineEVALUATE_TP_ENERG_DIFF, inlineEVALUATE_TP_ENERG_DIFF)
#define FORTNT_EVALUATE_TP_ENERG_DIFF FORTRAN_NAME( EVALUATE_TP_ENERG_DIFF, evaluate_tp_energ_diff)

inline void FORTRAN_NAME(inlineEVALUATE_TP_ENERG_DIFF, inlineEVALUATE_TP_ENERG_DIFF)(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m) )
{
 CH_TIMELEAF("FORT_EVALUATE_TP_ENERG_DIFF");
 FORTRAN_NAME( EVALUATE_TP_ENERG_DIFF ,evaluate_tp_energ_diff )(
      CHFt_FRA(ClsFlux)
      ,CHFt_CONST_FRA1(f)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(m) );
}
#endif  // GUARDEVALUATE_TP_ENERG_DIFF 

#ifndef GUARDEVALUATE_TP_ENERG_DIFF_SPECIES 
#define GUARDEVALUATE_TP_ENERG_DIFF_SPECIES 
// Prototype for Fortran procedure evaluate_tp_energ_diff_species ...
//
void FORTRAN_NAME( EVALUATE_TP_ENERG_DIFF_SPECIES ,evaluate_tp_energ_diff_species )(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m)
      ,CHFp_CONST_REAL(mb) );

#define FORT_EVALUATE_TP_ENERG_DIFF_SPECIES FORTRAN_NAME( inlineEVALUATE_TP_ENERG_DIFF_SPECIES, inlineEVALUATE_TP_ENERG_DIFF_SPECIES)
#define FORTNT_EVALUATE_TP_ENERG_DIFF_SPECIES FORTRAN_NAME( EVALUATE_TP_ENERG_DIFF_SPECIES, evaluate_tp_energ_diff_species)

inline void FORTRAN_NAME(inlineEVALUATE_TP_ENERG_DIFF_SPECIES, inlineEVALUATE_TP_ENERG_DIFF_SPECIES)(
      CHFp_FRA(ClsFlux)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(m)
      ,CHFp_CONST_REAL(mb) )
{
 CH_TIMELEAF("FORT_EVALUATE_TP_ENERG_DIFF_SPECIES");
 FORTRAN_NAME( EVALUATE_TP_ENERG_DIFF_SPECIES ,evaluate_tp_energ_diff_species )(
      CHFt_FRA(ClsFlux)
      ,CHFt_CONST_FRA1(f)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(m)
      ,CHFt_CONST_REAL(mb) );
}
#endif  // GUARDEVALUATE_TP_ENERG_DIFF_SPECIES 

#ifndef GUARDEVALUATE_COLL_FLUX_COMBINED 
#define GUARDEVALUATE_COLL_FLUX_COMBINED 
// Prototype for Fortran procedure evaluate_coll_flux_combined ...
//
void FORTRAN_NAME( EVALUATE_COLL_FLUX_COMBINED ,evaluate_coll_flux_combined )(
      CHFp_FRA(fluxRHS)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(fluxFace)
      ,CHFp_CONST_FRA(fluxCell)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) );

#define FORT_EVALUATE_COLL_FLUX_COMBINED FORTRAN_NAME( inlineEVALUATE_COLL_FLUX_COMBINED, inlineEVALUATE_COLL_FLUX_COMBINED)
#define FORTNT_EVALUATE_COLL_FLUX_COMBINED FORTRAN_NAME( EVALUATE_COLL_FLUX_COMBINED, evaluate_coll_flux_combined)

inline void FORTRAN_NAME(inlineEVALUATE_COLL_FLUX_COMBINED, inlineEVALUATE_COLL_FLUX_COMBINED)(
      CHFp_FRA(fluxRHS)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(fluxFace)
      ,CHFp_CONST_FRA(fluxCell)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu) )
{
 CH_TIMELEAF("FORT_EVALUATE_COLL_FLUX_COMBINED");
 FORTRAN_NAME( EVALUATE_COLL_FLUX_COMBINED ,evaluate_coll_flux_combined )(
      CHFt_FRA(fluxRHS)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(fluxFace)
      ,CHFt_CONST_FRA(fluxCell)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu) );
}
#endif  // GUARDEVALUATE_COLL_FLUX_COMBINED 

#ifndef GUARDEVALUATE_FULL_ER_FLUX 
#define GUARDEVALUATE_FULL_ER_FLUX 
// Prototype for Fortran procedure evaluate_full_er_flux ...
//
void FORTRAN_NAME( EVALUATE_FULL_ER_FLUX ,evaluate_full_er_flux )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(ERest)
      ,CHFp_CONST_FRA1(ENorm)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(m) );

#define FORT_EVALUATE_FULL_ER_FLUX FORTRAN_NAME( inlineEVALUATE_FULL_ER_FLUX, inlineEVALUATE_FULL_ER_FLUX)
#define FORTNT_EVALUATE_FULL_ER_FLUX FORTRAN_NAME( EVALUATE_FULL_ER_FLUX, evaluate_full_er_flux)

inline void FORTRAN_NAME(inlineEVALUATE_FULL_ER_FLUX, inlineEVALUATE_FULL_ER_FLUX)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(ERest)
      ,CHFp_CONST_FRA1(ENorm)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(m) )
{
 CH_TIMELEAF("FORT_EVALUATE_FULL_ER_FLUX");
 FORTRAN_NAME( EVALUATE_FULL_ER_FLUX ,evaluate_full_er_flux )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(box)
      ,CHFt_FRA(flux)
      ,CHFt_CONST_FRA1(ERest)
      ,CHFt_CONST_FRA1(ENorm)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu)
      ,CHFt_CONST_REAL(m) );
}
#endif  // GUARDEVALUATE_FULL_ER_FLUX 

#ifndef GUARDEVALUATE_NORM_ER_FLUX 
#define GUARDEVALUATE_NORM_ER_FLUX 
// Prototype for Fortran procedure evaluate_norm_er_flux ...
//
void FORTRAN_NAME( EVALUATE_NORM_ER_FLUX ,evaluate_norm_er_flux )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(m) );

#define FORT_EVALUATE_NORM_ER_FLUX FORTRAN_NAME( inlineEVALUATE_NORM_ER_FLUX, inlineEVALUATE_NORM_ER_FLUX)
#define FORTNT_EVALUATE_NORM_ER_FLUX FORTRAN_NAME( EVALUATE_NORM_ER_FLUX, evaluate_norm_er_flux)

inline void FORTRAN_NAME(inlineEVALUATE_NORM_ER_FLUX, inlineEVALUATE_NORM_ER_FLUX)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(m) )
{
 CH_TIMELEAF("FORT_EVALUATE_NORM_ER_FLUX");
 FORTRAN_NAME( EVALUATE_NORM_ER_FLUX ,evaluate_norm_er_flux )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(box)
      ,CHFt_FRA(flux)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu)
      ,CHFt_CONST_REAL(m) );
}
#endif  // GUARDEVALUATE_NORM_ER_FLUX 

#ifndef GUARDCOMPUTE_SC_CLS_FREQ 
#define GUARDCOMPUTE_SC_CLS_FREQ 
// Prototype for Fortran procedure compute_sc_cls_freq ...
//
void FORTRAN_NAME( COMPUTE_SC_CLS_FREQ ,compute_sc_cls_freq )(
      CHFp_BOX(box)
      ,CHFp_FRA1(cls_freq)
      ,CHFp_CONST_FRA1(n)
      ,CHFp_CONST_FRA1(T) );

#define FORT_COMPUTE_SC_CLS_FREQ FORTRAN_NAME( inlineCOMPUTE_SC_CLS_FREQ, inlineCOMPUTE_SC_CLS_FREQ)
#define FORTNT_COMPUTE_SC_CLS_FREQ FORTRAN_NAME( COMPUTE_SC_CLS_FREQ, compute_sc_cls_freq)

inline void FORTRAN_NAME(inlineCOMPUTE_SC_CLS_FREQ, inlineCOMPUTE_SC_CLS_FREQ)(
      CHFp_BOX(box)
      ,CHFp_FRA1(cls_freq)
      ,CHFp_CONST_FRA1(n)
      ,CHFp_CONST_FRA1(T) )
{
 CH_TIMELEAF("FORT_COMPUTE_SC_CLS_FREQ");
 FORTRAN_NAME( COMPUTE_SC_CLS_FREQ ,compute_sc_cls_freq )(
      CHFt_BOX(box)
      ,CHFt_FRA1(cls_freq)
      ,CHFt_CONST_FRA1(n)
      ,CHFt_CONST_FRA1(T) );
}
#endif  // GUARDCOMPUTE_SC_CLS_FREQ 

#ifndef GUARDEVAL_CONSDRAGDIFF_FLUX 
#define GUARDEVAL_CONSDRAGDIFF_FLUX 
// Prototype for Fortran procedure eval_consdragdiff_flux ...
//
void FORTRAN_NAME( EVAL_CONSDRAGDIFF_FLUX ,eval_consdragdiff_flux )(
      CHFp_BOX(box)
      ,CHFp_FRA(fluxes)
      ,CHFp_CONST_FRA1(fBJ)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(Nvp)
      ,CHFp_CONST_INT(Nmu) );

#define FORT_EVAL_CONSDRAGDIFF_FLUX FORTRAN_NAME( inlineEVAL_CONSDRAGDIFF_FLUX, inlineEVAL_CONSDRAGDIFF_FLUX)
#define FORTNT_EVAL_CONSDRAGDIFF_FLUX FORTRAN_NAME( EVAL_CONSDRAGDIFF_FLUX, eval_consdragdiff_flux)

inline void FORTRAN_NAME(inlineEVAL_CONSDRAGDIFF_FLUX, inlineEVAL_CONSDRAGDIFF_FLUX)(
      CHFp_BOX(box)
      ,CHFp_FRA(fluxes)
      ,CHFp_CONST_FRA1(fBJ)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(Nvp)
      ,CHFp_CONST_INT(Nmu) )
{
 CH_TIMELEAF("FORT_EVAL_CONSDRAGDIFF_FLUX");
 FORTRAN_NAME( EVAL_CONSDRAGDIFF_FLUX ,eval_consdragdiff_flux )(
      CHFt_BOX(box)
      ,CHFt_FRA(fluxes)
      ,CHFt_CONST_FRA1(fBJ)
      ,CHFt_CONST_FRA1(B)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(Nvp)
      ,CHFt_CONST_INT(Nmu) );
}
#endif  // GUARDEVAL_CONSDRAGDIFF_FLUX 

#ifndef GUARDEVAL_CONS_UPAR_TEMP 
#define GUARDEVAL_CONS_UPAR_TEMP 
// Prototype for Fortran procedure eval_cons_upar_temp ...
//
void FORTRAN_NAME( EVAL_CONS_UPAR_TEMP ,eval_cons_upar_temp )(
      CHFp_BOX(box)
      ,CHFp_FRA1(Upar)
      ,CHFp_FRA1(Temp)
      ,CHFp_CONST_FRA(vmoms)
      ,CHFp_CONST_FRA(pmoms) );

#define FORT_EVAL_CONS_UPAR_TEMP FORTRAN_NAME( inlineEVAL_CONS_UPAR_TEMP, inlineEVAL_CONS_UPAR_TEMP)
#define FORTNT_EVAL_CONS_UPAR_TEMP FORTRAN_NAME( EVAL_CONS_UPAR_TEMP, eval_cons_upar_temp)

inline void FORTRAN_NAME(inlineEVAL_CONS_UPAR_TEMP, inlineEVAL_CONS_UPAR_TEMP)(
      CHFp_BOX(box)
      ,CHFp_FRA1(Upar)
      ,CHFp_FRA1(Temp)
      ,CHFp_CONST_FRA(vmoms)
      ,CHFp_CONST_FRA(pmoms) )
{
 CH_TIMELEAF("FORT_EVAL_CONS_UPAR_TEMP");
 FORTRAN_NAME( EVAL_CONS_UPAR_TEMP ,eval_cons_upar_temp )(
      CHFt_BOX(box)
      ,CHFt_FRA1(Upar)
      ,CHFt_FRA1(Temp)
      ,CHFt_CONST_FRA(vmoms)
      ,CHFt_CONST_FRA(pmoms) );
}
#endif  // GUARDEVAL_CONS_UPAR_TEMP 

#ifndef GUARDEVAL_CONSDRAGDIFF_RHS 
#define GUARDEVAL_CONSDRAGDIFF_RHS 
// Prototype for Fortran procedure eval_consdragdiff_rhs ...
//
void FORTRAN_NAME( EVAL_CONSDRAGDIFF_RHS ,eval_consdragdiff_rhs )(
      CHFp_BOX(box)
      ,CHFp_FRA1(rhs_cls)
      ,CHFp_CONST_REAL(nu)
      ,CHFp_CONST_FRA1(Upar)
      ,CHFp_CONST_FRA1(Temp)
      ,CHFp_CONST_FRA(Jpsi) );

#define FORT_EVAL_CONSDRAGDIFF_RHS FORTRAN_NAME( inlineEVAL_CONSDRAGDIFF_RHS, inlineEVAL_CONSDRAGDIFF_RHS)
#define FORTNT_EVAL_CONSDRAGDIFF_RHS FORTRAN_NAME( EVAL_CONSDRAGDIFF_RHS, eval_consdragdiff_rhs)

inline void FORTRAN_NAME(inlineEVAL_CONSDRAGDIFF_RHS, inlineEVAL_CONSDRAGDIFF_RHS)(
      CHFp_BOX(box)
      ,CHFp_FRA1(rhs_cls)
      ,CHFp_CONST_REAL(nu)
      ,CHFp_CONST_FRA1(Upar)
      ,CHFp_CONST_FRA1(Temp)
      ,CHFp_CONST_FRA(Jpsi) )
{
 CH_TIMELEAF("FORT_EVAL_CONSDRAGDIFF_RHS");
 FORTRAN_NAME( EVAL_CONSDRAGDIFF_RHS ,eval_consdragdiff_rhs )(
      CHFt_BOX(box)
      ,CHFt_FRA1(rhs_cls)
      ,CHFt_CONST_REAL(nu)
      ,CHFt_CONST_FRA1(Upar)
      ,CHFt_CONST_FRA1(Temp)
      ,CHFt_CONST_FRA(Jpsi) );
}
#endif  // GUARDEVAL_CONSDRAGDIFF_RHS 

}

#endif
