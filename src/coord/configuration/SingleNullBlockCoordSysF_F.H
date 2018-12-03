#ifndef _SINGLENULLBLOCKCOORDSYSF_F_H_
#define _SINGLENULLBLOCKCOORDSYSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGET_FIELD_FROM_RBPOL_RBTOR 
#define GUARDGET_FIELD_FROM_RBPOL_RBTOR 
// Prototype for Fortran procedure get_field_from_rbpol_rbtor ...
//
void FORTRAN_NAME( GET_FIELD_FROM_RBPOL_RBTOR ,get_field_from_rbpol_rbtor )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_CONST_FRA(RBpol)
      ,CHFp_CONST_FRA1(R)
      ,CHFp_FRA(B) );

#define FORT_GET_FIELD_FROM_RBPOL_RBTOR FORTRAN_NAME( inlineGET_FIELD_FROM_RBPOL_RBTOR, inlineGET_FIELD_FROM_RBPOL_RBTOR)
#define FORTNT_GET_FIELD_FROM_RBPOL_RBTOR FORTRAN_NAME( GET_FIELD_FROM_RBPOL_RBTOR, get_field_from_rbpol_rbtor)

inline void FORTRAN_NAME(inlineGET_FIELD_FROM_RBPOL_RBTOR, inlineGET_FIELD_FROM_RBPOL_RBTOR)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_CONST_FRA(RBpol)
      ,CHFp_CONST_FRA1(R)
      ,CHFp_FRA(B) )
{
 CH_TIMELEAF("FORT_GET_FIELD_FROM_RBPOL_RBTOR");
 FORTRAN_NAME( GET_FIELD_FROM_RBPOL_RBTOR ,get_field_from_rbpol_rbtor )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_REAL(RBtor)
      ,CHFt_CONST_FRA(RBpol)
      ,CHFt_CONST_FRA1(R)
      ,CHFt_FRA(B) );
}
#endif  // GUARDGET_FIELD_FROM_RBPOL_RBTOR 

#ifndef GUARDGET_FIELD_MAGNITUDE_AND_UNITVECTOR 
#define GUARDGET_FIELD_MAGNITUDE_AND_UNITVECTOR 
// Prototype for Fortran procedure get_field_magnitude_and_unitvector ...
//
void FORTRAN_NAME( GET_FIELD_MAGNITUDE_AND_UNITVECTOR ,get_field_magnitude_and_unitvector )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(Bvec)
      ,CHFp_FRA1(Bmag)
      ,CHFp_FRA(bunit) );

#define FORT_GET_FIELD_MAGNITUDE_AND_UNITVECTOR FORTRAN_NAME( inlineGET_FIELD_MAGNITUDE_AND_UNITVECTOR, inlineGET_FIELD_MAGNITUDE_AND_UNITVECTOR)
#define FORTNT_GET_FIELD_MAGNITUDE_AND_UNITVECTOR FORTRAN_NAME( GET_FIELD_MAGNITUDE_AND_UNITVECTOR, get_field_magnitude_and_unitvector)

inline void FORTRAN_NAME(inlineGET_FIELD_MAGNITUDE_AND_UNITVECTOR, inlineGET_FIELD_MAGNITUDE_AND_UNITVECTOR)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(Bvec)
      ,CHFp_FRA1(Bmag)
      ,CHFp_FRA(bunit) )
{
 CH_TIMELEAF("FORT_GET_FIELD_MAGNITUDE_AND_UNITVECTOR");
 FORTRAN_NAME( GET_FIELD_MAGNITUDE_AND_UNITVECTOR ,get_field_magnitude_and_unitvector )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(Bvec)
      ,CHFt_FRA1(Bmag)
      ,CHFt_FRA(bunit) );
}
#endif  // GUARDGET_FIELD_MAGNITUDE_AND_UNITVECTOR 

#ifndef GUARDGET_FIELD_DERIVATIVE_DATA 
#define GUARDGET_FIELD_DERIVATIVE_DATA 
// Prototype for Fortran procedure get_field_derivative_data ...
//
void FORTRAN_NAME( GET_FIELD_DERIVATIVE_DATA ,get_field_derivative_data )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_INT(axisymmetric)
      ,CHFp_CONST_FRA1(RB)
      ,CHFp_CONST_FRA1(dRBdr)
      ,CHFp_CONST_FRA1(dRBdt)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA1(dbunitRdr)
      ,CHFp_CONST_FRA1(dbunitRdt)
      ,CHFp_CONST_FRA1(dbunitphidr)
      ,CHFp_CONST_FRA1(dbunitphidt)
      ,CHFp_CONST_FRA1(dbunitZdr)
      ,CHFp_CONST_FRA1(dbunitZdt)
      ,CHFp_CONST_FRA1(R)
      ,CHFp_CONST_FRA1(Rr)
      ,CHFp_CONST_FRA1(Rt)
      ,CHFp_CONST_FRA1(Zr)
      ,CHFp_CONST_FRA1(Zt)
      ,CHFp_FRA(gradb)
      ,CHFp_FRA(curlb)
      ,CHFp_FRA1(bdotcurlb) );

#define FORT_GET_FIELD_DERIVATIVE_DATA FORTRAN_NAME( inlineGET_FIELD_DERIVATIVE_DATA, inlineGET_FIELD_DERIVATIVE_DATA)
#define FORTNT_GET_FIELD_DERIVATIVE_DATA FORTRAN_NAME( GET_FIELD_DERIVATIVE_DATA, get_field_derivative_data)

inline void FORTRAN_NAME(inlineGET_FIELD_DERIVATIVE_DATA, inlineGET_FIELD_DERIVATIVE_DATA)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_INT(axisymmetric)
      ,CHFp_CONST_FRA1(RB)
      ,CHFp_CONST_FRA1(dRBdr)
      ,CHFp_CONST_FRA1(dRBdt)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA1(dbunitRdr)
      ,CHFp_CONST_FRA1(dbunitRdt)
      ,CHFp_CONST_FRA1(dbunitphidr)
      ,CHFp_CONST_FRA1(dbunitphidt)
      ,CHFp_CONST_FRA1(dbunitZdr)
      ,CHFp_CONST_FRA1(dbunitZdt)
      ,CHFp_CONST_FRA1(R)
      ,CHFp_CONST_FRA1(Rr)
      ,CHFp_CONST_FRA1(Rt)
      ,CHFp_CONST_FRA1(Zr)
      ,CHFp_CONST_FRA1(Zt)
      ,CHFp_FRA(gradb)
      ,CHFp_FRA(curlb)
      ,CHFp_FRA1(bdotcurlb) )
{
 CH_TIMELEAF("FORT_GET_FIELD_DERIVATIVE_DATA");
 FORTRAN_NAME( GET_FIELD_DERIVATIVE_DATA ,get_field_derivative_data )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_INT(axisymmetric)
      ,CHFt_CONST_FRA1(RB)
      ,CHFt_CONST_FRA1(dRBdr)
      ,CHFt_CONST_FRA1(dRBdt)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_CONST_FRA1(dbunitRdr)
      ,CHFt_CONST_FRA1(dbunitRdt)
      ,CHFt_CONST_FRA1(dbunitphidr)
      ,CHFt_CONST_FRA1(dbunitphidt)
      ,CHFt_CONST_FRA1(dbunitZdr)
      ,CHFt_CONST_FRA1(dbunitZdt)
      ,CHFt_CONST_FRA1(R)
      ,CHFt_CONST_FRA1(Rr)
      ,CHFt_CONST_FRA1(Rt)
      ,CHFt_CONST_FRA1(Zr)
      ,CHFt_CONST_FRA1(Zt)
      ,CHFt_FRA(gradb)
      ,CHFt_FRA(curlb)
      ,CHFt_FRA1(bdotcurlb) );
}
#endif  // GUARDGET_FIELD_DERIVATIVE_DATA 

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
