#ifndef _FPLIMITERSF_F_H_
#define _FPLIMITERSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRECONSTRUCT_DFN_FACE_VPAR 
#define GUARDRECONSTRUCT_DFN_FACE_VPAR 
// Prototype for Fortran procedure reconstruct_dfn_face_vpar ...
//
void FORTRAN_NAME( RECONSTRUCT_DFN_FACE_VPAR ,reconstruct_dfn_face_vpar )(
      CHFp_FRA1(f)
      ,CHFp_CONST_FRA1(a)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_INT(flag)
      ,CHFp_CONST_INT(Nvpar) );

#define FORT_RECONSTRUCT_DFN_FACE_VPAR FORTRAN_NAME( inlineRECONSTRUCT_DFN_FACE_VPAR, inlineRECONSTRUCT_DFN_FACE_VPAR)
#define FORTNT_RECONSTRUCT_DFN_FACE_VPAR FORTRAN_NAME( RECONSTRUCT_DFN_FACE_VPAR, reconstruct_dfn_face_vpar)

inline void FORTRAN_NAME(inlineRECONSTRUCT_DFN_FACE_VPAR, inlineRECONSTRUCT_DFN_FACE_VPAR)(
      CHFp_FRA1(f)
      ,CHFp_CONST_FRA1(a)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_INT(flag)
      ,CHFp_CONST_INT(Nvpar) )
{
 CH_TIMELEAF("FORT_RECONSTRUCT_DFN_FACE_VPAR");
 FORTRAN_NAME( RECONSTRUCT_DFN_FACE_VPAR ,reconstruct_dfn_face_vpar )(
      CHFt_FRA1(f)
      ,CHFt_CONST_FRA1(a)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_INT(flag)
      ,CHFt_CONST_INT(Nvpar) );
}
#endif  // GUARDRECONSTRUCT_DFN_FACE_VPAR 

#ifndef GUARDRECONSTRUCT_DFN_FACE_MU 
#define GUARDRECONSTRUCT_DFN_FACE_MU 
// Prototype for Fortran procedure reconstruct_dfn_face_mu ...
//
void FORTRAN_NAME( RECONSTRUCT_DFN_FACE_MU ,reconstruct_dfn_face_mu )(
      CHFp_FRA1(f)
      ,CHFp_CONST_FRA1(a)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_INT(flag)
      ,CHFp_CONST_INT(Nmu) );

#define FORT_RECONSTRUCT_DFN_FACE_MU FORTRAN_NAME( inlineRECONSTRUCT_DFN_FACE_MU, inlineRECONSTRUCT_DFN_FACE_MU)
#define FORTNT_RECONSTRUCT_DFN_FACE_MU FORTRAN_NAME( RECONSTRUCT_DFN_FACE_MU, reconstruct_dfn_face_mu)

inline void FORTRAN_NAME(inlineRECONSTRUCT_DFN_FACE_MU, inlineRECONSTRUCT_DFN_FACE_MU)(
      CHFp_FRA1(f)
      ,CHFp_CONST_FRA1(a)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_INT(flag)
      ,CHFp_CONST_INT(Nmu) )
{
 CH_TIMELEAF("FORT_RECONSTRUCT_DFN_FACE_MU");
 FORTRAN_NAME( RECONSTRUCT_DFN_FACE_MU ,reconstruct_dfn_face_mu )(
      CHFt_FRA1(f)
      ,CHFt_CONST_FRA1(a)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_INT(flag)
      ,CHFt_CONST_INT(Nmu) );
}
#endif  // GUARDRECONSTRUCT_DFN_FACE_MU 

}

#endif
