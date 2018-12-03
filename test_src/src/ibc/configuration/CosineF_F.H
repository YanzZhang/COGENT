#ifndef _COSINEF_F_H_
#define _COSINEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSET_COSINE 
#define GUARDSET_COSINE 
// Prototype for Fortran procedure set_cosine ...
//
void FORTRAN_NAME( SET_COSINE ,set_cosine )(
      CHFp_FRA(phi)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(shift)
      ,CHFp_CONST_REAL(constant)
      ,CHFp_CONST_REAL(amplitude)
      ,CHFp_CONST_REALVECT(mode)
      ,CHFp_CONST_REALVECT(phase) );

#define FORT_SET_COSINE FORTRAN_NAME( inlineSET_COSINE, inlineSET_COSINE)
#define FORTNT_SET_COSINE FORTRAN_NAME( SET_COSINE, set_cosine)

inline void FORTRAN_NAME(inlineSET_COSINE, inlineSET_COSINE)(
      CHFp_FRA(phi)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(shift)
      ,CHFp_CONST_REAL(constant)
      ,CHFp_CONST_REAL(amplitude)
      ,CHFp_CONST_REALVECT(mode)
      ,CHFp_CONST_REALVECT(phase) )
{
 CH_TIMELEAF("FORT_SET_COSINE");
 FORTRAN_NAME( SET_COSINE ,set_cosine )(
      CHFt_FRA(phi)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(x)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(shift)
      ,CHFt_CONST_REAL(constant)
      ,CHFt_CONST_REAL(amplitude)
      ,CHFt_CONST_REALVECT(mode)
      ,CHFt_CONST_REALVECT(phase) );
}
#endif  // GUARDSET_COSINE 

}

#endif
