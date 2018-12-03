#ifndef _PHASEBLOCKCOORDSYSF_F_H_
#define _PHASEBLOCKCOORDSYSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDPHASE_BLOCK_PROJECT 
#define GUARDPHASE_BLOCK_PROJECT 
// Prototype for Fortran procedure phase_block_project ...
//
void FORTRAN_NAME( PHASE_BLOCK_PROJECT ,phase_block_project )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(vec_src)
      ,CHFp_FRA(vec_dst) );

#define FORT_PHASE_BLOCK_PROJECT FORTRAN_NAME( inlinePHASE_BLOCK_PROJECT, inlinePHASE_BLOCK_PROJECT)
#define FORTNT_PHASE_BLOCK_PROJECT FORTRAN_NAME( PHASE_BLOCK_PROJECT, phase_block_project)

inline void FORTRAN_NAME(inlinePHASE_BLOCK_PROJECT, inlinePHASE_BLOCK_PROJECT)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(vec_src)
      ,CHFp_FRA(vec_dst) )
{
 CH_TIMELEAF("FORT_PHASE_BLOCK_PROJECT");
 FORTRAN_NAME( PHASE_BLOCK_PROJECT ,phase_block_project )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(vec_src)
      ,CHFt_FRA(vec_dst) );
}
#endif  // GUARDPHASE_BLOCK_PROJECT 

#ifndef GUARDPHASE_BLOCK_PSITHETA_PROJECTIONS 
#define GUARDPHASE_BLOCK_PSITHETA_PROJECTIONS 
// Prototype for Fortran procedure phase_block_psitheta_projections ...
//
void FORTRAN_NAME( PHASE_BLOCK_PSITHETA_PROJECTIONS ,phase_block_psitheta_projections )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(vec_psi)
      ,CHFp_CONST_FRA(vec_theta)
      ,CHFp_FRA(vec_dst) );

#define FORT_PHASE_BLOCK_PSITHETA_PROJECTIONS FORTRAN_NAME( inlinePHASE_BLOCK_PSITHETA_PROJECTIONS, inlinePHASE_BLOCK_PSITHETA_PROJECTIONS)
#define FORTNT_PHASE_BLOCK_PSITHETA_PROJECTIONS FORTRAN_NAME( PHASE_BLOCK_PSITHETA_PROJECTIONS, phase_block_psitheta_projections)

inline void FORTRAN_NAME(inlinePHASE_BLOCK_PSITHETA_PROJECTIONS, inlinePHASE_BLOCK_PSITHETA_PROJECTIONS)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(vec_psi)
      ,CHFp_CONST_FRA(vec_theta)
      ,CHFp_FRA(vec_dst) )
{
 CH_TIMELEAF("FORT_PHASE_BLOCK_PSITHETA_PROJECTIONS");
 FORTRAN_NAME( PHASE_BLOCK_PSITHETA_PROJECTIONS ,phase_block_psitheta_projections )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(vec_psi)
      ,CHFt_CONST_FRA(vec_theta)
      ,CHFt_FRA(vec_dst) );
}
#endif  // GUARDPHASE_BLOCK_PSITHETA_PROJECTIONS 

#ifndef GUARDPHASE_BLOCK_GRADF_FACTOR 
#define GUARDPHASE_BLOCK_GRADF_FACTOR 
// Prototype for Fortran procedure phase_block_gradf_factor ...
//
void FORTRAN_NAME( PHASE_BLOCK_GRADF_FACTOR ,phase_block_gradf_factor )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA1(vec_src)
      ,CHFp_FRA1(vec_dst) );

#define FORT_PHASE_BLOCK_GRADF_FACTOR FORTRAN_NAME( inlinePHASE_BLOCK_GRADF_FACTOR, inlinePHASE_BLOCK_GRADF_FACTOR)
#define FORTNT_PHASE_BLOCK_GRADF_FACTOR FORTRAN_NAME( PHASE_BLOCK_GRADF_FACTOR, phase_block_gradf_factor)

inline void FORTRAN_NAME(inlinePHASE_BLOCK_GRADF_FACTOR, inlinePHASE_BLOCK_GRADF_FACTOR)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA1(vec_src)
      ,CHFp_FRA1(vec_dst) )
{
 CH_TIMELEAF("FORT_PHASE_BLOCK_GRADF_FACTOR");
 FORTRAN_NAME( PHASE_BLOCK_GRADF_FACTOR ,phase_block_gradf_factor )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA1(vec_src)
      ,CHFt_FRA1(vec_dst) );
}
#endif  // GUARDPHASE_BLOCK_GRADF_FACTOR 

}

#endif
