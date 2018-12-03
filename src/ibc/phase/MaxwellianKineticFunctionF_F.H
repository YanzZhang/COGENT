#ifndef _MAXWELLIANKINETICFUNCTIONF_F_H_
#define _MAXWELLIANKINETICFUNCTIONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSET_MAXWELL4D 
#define GUARDSET_MAXWELL4D 
// Prototype for Fortran procedure set_maxwell4d ...
//
void FORTRAN_NAME( SET_MAXWELL4D ,set_maxwell4d )(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coords)
      ,CHFp_CONST_FRA1(dens)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_CONST_FRA1(vshift)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_REAL(mass) );

#define FORT_SET_MAXWELL4D FORTRAN_NAME( inlineSET_MAXWELL4D, inlineSET_MAXWELL4D)
#define FORTNT_SET_MAXWELL4D FORTRAN_NAME( SET_MAXWELL4D, set_maxwell4d)

inline void FORTRAN_NAME(inlineSET_MAXWELL4D, inlineSET_MAXWELL4D)(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coords)
      ,CHFp_CONST_FRA1(dens)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_CONST_FRA1(vshift)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_REAL(mass) )
{
 CH_TIMELEAF("FORT_SET_MAXWELL4D");
 FORTRAN_NAME( SET_MAXWELL4D ,set_maxwell4d )(
      CHFt_FRA(f)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(coords)
      ,CHFt_CONST_FRA1(dens)
      ,CHFt_CONST_FRA1(temp)
      ,CHFt_CONST_FRA1(vshift)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_CONST_REAL(mass) );
}
#endif  // GUARDSET_MAXWELL4D 

}

#endif
