#ifndef _GKBOLTZMANNRELATIONF_F_H_
#define _GKBOLTZMANNRELATIONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_ELECTRON_DENSITY 
#define GUARDCOMPUTE_ELECTRON_DENSITY 
// Prototype for Fortran procedure compute_electron_density ...
//
void FORTRAN_NAME( COMPUTE_ELECTRON_DENSITY ,compute_electron_density )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(te)
      ,CHFp_FRA1(ne) );

#define FORT_COMPUTE_ELECTRON_DENSITY FORTRAN_NAME( inlineCOMPUTE_ELECTRON_DENSITY, inlineCOMPUTE_ELECTRON_DENSITY)
#define FORTNT_COMPUTE_ELECTRON_DENSITY FORTRAN_NAME( COMPUTE_ELECTRON_DENSITY, compute_electron_density)

inline void FORTRAN_NAME(inlineCOMPUTE_ELECTRON_DENSITY, inlineCOMPUTE_ELECTRON_DENSITY)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(te)
      ,CHFp_FRA1(ne) )
{
 CH_TIMELEAF("FORT_COMPUTE_ELECTRON_DENSITY");
 FORTRAN_NAME( COMPUTE_ELECTRON_DENSITY ,compute_electron_density )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(te)
      ,CHFt_FRA1(ne) );
}
#endif  // GUARDCOMPUTE_ELECTRON_DENSITY 

}

#endif
