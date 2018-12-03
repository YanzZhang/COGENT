#ifndef _CODIMBCF_F_H_
#define _CODIMBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_LIMITS 
#define GUARDCOMPUTE_LIMITS 
// Prototype for Fortran procedure compute_limits ...
//
void FORTRAN_NAME( COMPUTE_LIMITS ,compute_limits )(
      CHFp_CONST_INTVECT(ibeg)
      ,CHFp_CONST_INTVECT(iend)
      ,CHFp_CONST_INTVECT(istride)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0)
      ,CHFp_CONST_INT(icodim) );

#define FORT_COMPUTE_LIMITS FORTRAN_NAME( inlineCOMPUTE_LIMITS, inlineCOMPUTE_LIMITS)
#define FORTNT_COMPUTE_LIMITS FORTRAN_NAME( COMPUTE_LIMITS, compute_limits)

inline void FORTRAN_NAME(inlineCOMPUTE_LIMITS, inlineCOMPUTE_LIMITS)(
      CHFp_CONST_INTVECT(ibeg)
      ,CHFp_CONST_INTVECT(iend)
      ,CHFp_CONST_INTVECT(istride)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0)
      ,CHFp_CONST_INT(icodim) )
{
 CH_TIMELEAF("FORT_COMPUTE_LIMITS");
 FORTRAN_NAME( COMPUTE_LIMITS ,compute_limits )(
      CHFt_CONST_INTVECT(ibeg)
      ,CHFt_CONST_INTVECT(iend)
      ,CHFt_CONST_INTVECT(istride)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_I1D(idir,iidirhi0)
      ,CHFt_CONST_I1D(iside,iisidehi0)
      ,CHFt_CONST_INT(icodim) );
}
#endif  // GUARDCOMPUTE_LIMITS 

#ifndef GUARDFILL_CODIM_GHOST_CELLS 
#define GUARDFILL_CODIM_GHOST_CELLS 
// Prototype for Fortran procedure fill_codim_ghost_cells ...
//
void FORTRAN_NAME( FILL_CODIM_GHOST_CELLS ,fill_codim_ghost_cells )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0)
      ,CHFp_CONST_INT(icodim) );

#define FORT_FILL_CODIM_GHOST_CELLS FORTRAN_NAME( inlineFILL_CODIM_GHOST_CELLS, inlineFILL_CODIM_GHOST_CELLS)
#define FORTNT_FILL_CODIM_GHOST_CELLS FORTRAN_NAME( FILL_CODIM_GHOST_CELLS, fill_codim_ghost_cells)

inline void FORTRAN_NAME(inlineFILL_CODIM_GHOST_CELLS, inlineFILL_CODIM_GHOST_CELLS)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0)
      ,CHFp_CONST_INT(icodim) )
{
 CH_TIMELEAF("FORT_FILL_CODIM_GHOST_CELLS");
 FORTRAN_NAME( FILL_CODIM_GHOST_CELLS ,fill_codim_ghost_cells )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_I1D(idir,iidirhi0)
      ,CHFt_CONST_I1D(iside,iisidehi0)
      ,CHFt_CONST_INT(icodim) );
}
#endif  // GUARDFILL_CODIM_GHOST_CELLS 

#ifndef GUARDFILL_CODIM2_GHOST_CELLS 
#define GUARDFILL_CODIM2_GHOST_CELLS 
// Prototype for Fortran procedure fill_codim2_ghost_cells ...
//
void FORTRAN_NAME( FILL_CODIM2_GHOST_CELLS ,fill_codim2_ghost_cells )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0) );

#define FORT_FILL_CODIM2_GHOST_CELLS FORTRAN_NAME( inlineFILL_CODIM2_GHOST_CELLS, inlineFILL_CODIM2_GHOST_CELLS)
#define FORTNT_FILL_CODIM2_GHOST_CELLS FORTRAN_NAME( FILL_CODIM2_GHOST_CELLS, fill_codim2_ghost_cells)

inline void FORTRAN_NAME(inlineFILL_CODIM2_GHOST_CELLS, inlineFILL_CODIM2_GHOST_CELLS)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0) )
{
 CH_TIMELEAF("FORT_FILL_CODIM2_GHOST_CELLS");
 FORTRAN_NAME( FILL_CODIM2_GHOST_CELLS ,fill_codim2_ghost_cells )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_I1D(idir,iidirhi0)
      ,CHFt_CONST_I1D(iside,iisidehi0) );
}
#endif  // GUARDFILL_CODIM2_GHOST_CELLS 

#ifndef GUARDFILL_CODIM3_GHOST_CELLS 
#define GUARDFILL_CODIM3_GHOST_CELLS 
// Prototype for Fortran procedure fill_codim3_ghost_cells ...
//
void FORTRAN_NAME( FILL_CODIM3_GHOST_CELLS ,fill_codim3_ghost_cells )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0) );

#define FORT_FILL_CODIM3_GHOST_CELLS FORTRAN_NAME( inlineFILL_CODIM3_GHOST_CELLS, inlineFILL_CODIM3_GHOST_CELLS)
#define FORTNT_FILL_CODIM3_GHOST_CELLS FORTRAN_NAME( FILL_CODIM3_GHOST_CELLS, fill_codim3_ghost_cells)

inline void FORTRAN_NAME(inlineFILL_CODIM3_GHOST_CELLS, inlineFILL_CODIM3_GHOST_CELLS)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0) )
{
 CH_TIMELEAF("FORT_FILL_CODIM3_GHOST_CELLS");
 FORTRAN_NAME( FILL_CODIM3_GHOST_CELLS ,fill_codim3_ghost_cells )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_I1D(idir,iidirhi0)
      ,CHFt_CONST_I1D(iside,iisidehi0) );
}
#endif  // GUARDFILL_CODIM3_GHOST_CELLS 

#ifndef GUARDFILL_CODIM4_GHOST_CELLS 
#define GUARDFILL_CODIM4_GHOST_CELLS 
// Prototype for Fortran procedure fill_codim4_ghost_cells ...
//
void FORTRAN_NAME( FILL_CODIM4_GHOST_CELLS ,fill_codim4_ghost_cells )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0) );

#define FORT_FILL_CODIM4_GHOST_CELLS FORTRAN_NAME( inlineFILL_CODIM4_GHOST_CELLS, inlineFILL_CODIM4_GHOST_CELLS)
#define FORTNT_FILL_CODIM4_GHOST_CELLS FORTRAN_NAME( FILL_CODIM4_GHOST_CELLS, fill_codim4_ghost_cells)

inline void FORTRAN_NAME(inlineFILL_CODIM4_GHOST_CELLS, inlineFILL_CODIM4_GHOST_CELLS)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_I1D(idir,iidirhi0)
      ,CHFp_CONST_I1D(iside,iisidehi0) )
{
 CH_TIMELEAF("FORT_FILL_CODIM4_GHOST_CELLS");
 FORTRAN_NAME( FILL_CODIM4_GHOST_CELLS ,fill_codim4_ghost_cells )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_I1D(idir,iidirhi0)
      ,CHFt_CONST_I1D(iside,iisidehi0) );
}
#endif  // GUARDFILL_CODIM4_GHOST_CELLS 

}

#endif
