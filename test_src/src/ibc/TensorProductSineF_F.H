#ifndef _TENSORPRODUCTSINEF_F_H_
#define _TENSORPRODUCTSINEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSET_TENSORPRODUCTSINE 
#define GUARDSET_TENSORPRODUCTSINE 
// Prototype for Fortran procedure set_tensorproductsine ...
//
void FORTRAN_NAME( SET_TENSORPRODUCTSINE ,set_tensorproductsine )(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coord)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(amp)
      ,CHFp_CONST_REALVECT(mode) );

#define FORT_SET_TENSORPRODUCTSINE FORTRAN_NAME( inlineSET_TENSORPRODUCTSINE, inlineSET_TENSORPRODUCTSINE)
#define FORTNT_SET_TENSORPRODUCTSINE FORTRAN_NAME( SET_TENSORPRODUCTSINE, set_tensorproductsine)

inline void FORTRAN_NAME(inlineSET_TENSORPRODUCTSINE, inlineSET_TENSORPRODUCTSINE)(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coord)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(amp)
      ,CHFp_CONST_REALVECT(mode) )
{
 CH_TIMELEAF("FORT_SET_TENSORPRODUCTSINE");
 FORTRAN_NAME( SET_TENSORPRODUCTSINE ,set_tensorproductsine )(
      CHFt_FRA(f)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(coord)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_REAL(amp)
      ,CHFt_CONST_REALVECT(mode) );
}
#endif  // GUARDSET_TENSORPRODUCTSINE 

}

#endif
