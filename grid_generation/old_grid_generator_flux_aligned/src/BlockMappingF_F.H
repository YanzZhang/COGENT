#ifndef _BLOCKMAPPINGF_F_H_
#define _BLOCKMAPPINGF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDDCT_INTERP 
#define GUARDDCT_INTERP 
// Prototype for Fortran procedure dct_interp ...
//
void FORTRAN_NAME( DCT_INTERP ,dct_interp )(
      CHFp_CONST_FRA1(coef)
      ,CHFp_CONST_INT(expansion_order)
      ,CHFp_CONST_INT(deriv1)
      ,CHFp_CONST_INT(deriv2)
      ,CHFp_R1D(fac1,ifac1hi0)
      ,CHFp_R1D(fac2,ifac2hi0)
      ,CHFp_R1D(sinfac1,isinfac1hi0)
      ,CHFp_R1D(cosfac1,icosfac1hi0)
      ,CHFp_R1D(sinfac2,isinfac2hi0)
      ,CHFp_R1D(cosfac2,icosfac2hi0)
      ,CHFp_R1D(lambda,ilambdahi0)
      ,CHFp_REAL(value) );

#define FORT_DCT_INTERP FORTRAN_NAME( inlineDCT_INTERP, inlineDCT_INTERP)
#define FORTNT_DCT_INTERP FORTRAN_NAME( DCT_INTERP, dct_interp)

inline void FORTRAN_NAME(inlineDCT_INTERP, inlineDCT_INTERP)(
      CHFp_CONST_FRA1(coef)
      ,CHFp_CONST_INT(expansion_order)
      ,CHFp_CONST_INT(deriv1)
      ,CHFp_CONST_INT(deriv2)
      ,CHFp_R1D(fac1,ifac1hi0)
      ,CHFp_R1D(fac2,ifac2hi0)
      ,CHFp_R1D(sinfac1,isinfac1hi0)
      ,CHFp_R1D(cosfac1,icosfac1hi0)
      ,CHFp_R1D(sinfac2,isinfac2hi0)
      ,CHFp_R1D(cosfac2,icosfac2hi0)
      ,CHFp_R1D(lambda,ilambdahi0)
      ,CHFp_REAL(value) )
{
 CH_TIMELEAF("FORT_DCT_INTERP");
 FORTRAN_NAME( DCT_INTERP ,dct_interp )(
      CHFt_CONST_FRA1(coef)
      ,CHFt_CONST_INT(expansion_order)
      ,CHFt_CONST_INT(deriv1)
      ,CHFt_CONST_INT(deriv2)
      ,CHFt_R1D(fac1,ifac1hi0)
      ,CHFt_R1D(fac2,ifac2hi0)
      ,CHFt_R1D(sinfac1,isinfac1hi0)
      ,CHFt_R1D(cosfac1,icosfac1hi0)
      ,CHFt_R1D(sinfac2,isinfac2hi0)
      ,CHFt_R1D(cosfac2,icosfac2hi0)
      ,CHFt_R1D(lambda,ilambdahi0)
      ,CHFt_REAL(value) );
}
#endif  // GUARDDCT_INTERP 

}

#endif
