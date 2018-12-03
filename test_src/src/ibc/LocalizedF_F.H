#ifndef _LOCALIZEDF_F_H_
#define _LOCALIZEDF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSET_LOCALIZED 
#define GUARDSET_LOCALIZED 
// Prototype for Fortran procedure set_localized ...
//
void FORTRAN_NAME( SET_LOCALIZED ,set_localized )(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coord)
      ,CHFp_CONST_REAL(amp)
      ,CHFp_CONST_REALVECT(center)
      ,CHFp_CONST_REALVECT(width)
      ,CHFp_CONST_REAL(floor) );

#define FORT_SET_LOCALIZED FORTRAN_NAME( inlineSET_LOCALIZED, inlineSET_LOCALIZED)
#define FORTNT_SET_LOCALIZED FORTRAN_NAME( SET_LOCALIZED, set_localized)

inline void FORTRAN_NAME(inlineSET_LOCALIZED, inlineSET_LOCALIZED)(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coord)
      ,CHFp_CONST_REAL(amp)
      ,CHFp_CONST_REALVECT(center)
      ,CHFp_CONST_REALVECT(width)
      ,CHFp_CONST_REAL(floor) )
{
 CH_TIMELEAF("FORT_SET_LOCALIZED");
 FORTRAN_NAME( SET_LOCALIZED ,set_localized )(
      CHFt_FRA(f)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(coord)
      ,CHFt_CONST_REAL(amp)
      ,CHFt_CONST_REALVECT(center)
      ,CHFt_CONST_REALVECT(width)
      ,CHFt_CONST_REAL(floor) );
}
#endif  // GUARDSET_LOCALIZED 

}

#endif
