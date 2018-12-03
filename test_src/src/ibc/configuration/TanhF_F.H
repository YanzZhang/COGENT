#ifndef _TANHF_F_H_
#define _TANHF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSET_TANH 
#define GUARDSET_TANH 
// Prototype for Fortran procedure set_tanh ...
//
void FORTRAN_NAME( SET_TANH ,set_tanh )(
      CHFp_FRA(phi)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(f1)
      ,CHFp_CONST_REAL(x1)
      ,CHFp_CONST_REAL(f2)
      ,CHFp_CONST_REAL(x2)
      ,CHFp_CONST_REAL(x0)
      ,CHFp_CONST_REAL(width) );

#define FORT_SET_TANH FORTRAN_NAME( inlineSET_TANH, inlineSET_TANH)
#define FORTNT_SET_TANH FORTRAN_NAME( SET_TANH, set_tanh)

inline void FORTRAN_NAME(inlineSET_TANH, inlineSET_TANH)(
      CHFp_FRA(phi)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(f1)
      ,CHFp_CONST_REAL(x1)
      ,CHFp_CONST_REAL(f2)
      ,CHFp_CONST_REAL(x2)
      ,CHFp_CONST_REAL(x0)
      ,CHFp_CONST_REAL(width) )
{
 CH_TIMELEAF("FORT_SET_TANH");
 FORTRAN_NAME( SET_TANH ,set_tanh )(
      CHFt_FRA(phi)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(x)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(f1)
      ,CHFt_CONST_REAL(x1)
      ,CHFt_CONST_REAL(f2)
      ,CHFt_CONST_REAL(x2)
      ,CHFt_CONST_REAL(x0)
      ,CHFt_CONST_REAL(width) );
}
#endif  // GUARDSET_TANH 

}

#endif
