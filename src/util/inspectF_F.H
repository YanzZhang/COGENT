#ifndef _INSPECTF_F_H_
#define _INSPECTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINSPECT 
#define GUARDINSPECT 
// Prototype for Fortran procedure inspect ...
//
void FORTRAN_NAME( INSPECT ,inspect )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(data) );

#define FORT_INSPECT FORTRAN_NAME( inlineINSPECT, inlineINSPECT)
#define FORTNT_INSPECT FORTRAN_NAME( INSPECT, inspect)

inline void FORTRAN_NAME(inlineINSPECT, inlineINSPECT)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(data) )
{
 CH_TIMELEAF("FORT_INSPECT");
 FORTRAN_NAME( INSPECT ,inspect )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(data) );
}
#endif  // GUARDINSPECT 

}

#endif
