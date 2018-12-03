#ifndef _MILLERBLOCKCOORDSYSF_F_H_
#define _MILLERBLOCKCOORDSYSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGET_MILLER_FIELD_DATA 
#define GUARDGET_MILLER_FIELD_DATA 
// Prototype for Fortran procedure get_miller_field_data ...
//
void FORTRAN_NAME( GET_MILLER_FIELD_DATA ,get_miller_field_data )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(RBpol)
      ,CHFp_CONST_FRA1(dRBpoldt)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_CONST_FRA1(Rmaj)
      ,CHFp_CONST_FRA1(Rr)
      ,CHFp_CONST_FRA1(Rt)
      ,CHFp_CONST_FRA1(Zr)
      ,CHFp_CONST_FRA1(Zt)
      ,CHFp_CONST_FRA1(Rtt)
      ,CHFp_CONST_FRA1(Ztt)
      ,CHFp_FRA(b_pt)
      ,CHFp_FRA1(Bmag_pt)
      ,CHFp_FRA(bunit_pt)
      ,CHFp_FRA(gradb_pt)
      ,CHFp_FRA(curlb_pt)
      ,CHFp_FRA1(bdotcurlb_pt) );

#define FORT_GET_MILLER_FIELD_DATA FORTRAN_NAME( inlineGET_MILLER_FIELD_DATA, inlineGET_MILLER_FIELD_DATA)
#define FORTNT_GET_MILLER_FIELD_DATA FORTRAN_NAME( GET_MILLER_FIELD_DATA, get_miller_field_data)

inline void FORTRAN_NAME(inlineGET_MILLER_FIELD_DATA, inlineGET_MILLER_FIELD_DATA)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(RBpol)
      ,CHFp_CONST_FRA1(dRBpoldt)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_CONST_FRA1(Rmaj)
      ,CHFp_CONST_FRA1(Rr)
      ,CHFp_CONST_FRA1(Rt)
      ,CHFp_CONST_FRA1(Zr)
      ,CHFp_CONST_FRA1(Zt)
      ,CHFp_CONST_FRA1(Rtt)
      ,CHFp_CONST_FRA1(Ztt)
      ,CHFp_FRA(b_pt)
      ,CHFp_FRA1(Bmag_pt)
      ,CHFp_FRA(bunit_pt)
      ,CHFp_FRA(gradb_pt)
      ,CHFp_FRA(curlb_pt)
      ,CHFp_FRA1(bdotcurlb_pt) )
{
 CH_TIMELEAF("FORT_GET_MILLER_FIELD_DATA");
 FORTRAN_NAME( GET_MILLER_FIELD_DATA ,get_miller_field_data )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(RBpol)
      ,CHFt_CONST_FRA1(dRBpoldt)
      ,CHFt_CONST_REAL(RBtor)
      ,CHFt_CONST_FRA1(Rmaj)
      ,CHFt_CONST_FRA1(Rr)
      ,CHFt_CONST_FRA1(Rt)
      ,CHFt_CONST_FRA1(Zr)
      ,CHFt_CONST_FRA1(Zt)
      ,CHFt_CONST_FRA1(Rtt)
      ,CHFt_CONST_FRA1(Ztt)
      ,CHFt_FRA(b_pt)
      ,CHFt_FRA1(Bmag_pt)
      ,CHFt_FRA(bunit_pt)
      ,CHFt_FRA(gradb_pt)
      ,CHFt_FRA(curlb_pt)
      ,CHFt_FRA1(bdotcurlb_pt) );
}
#endif  // GUARDGET_MILLER_FIELD_DATA 

#ifndef GUARDGET_RBPOL_MILLER 
#define GUARDGET_RBPOL_MILLER 
// Prototype for Fortran procedure get_rbpol_miller ...
//
void FORTRAN_NAME( GET_RBPOL_MILLER ,get_rbpol_miller )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(kappa)
      ,CHFp_CONST_REAL(dpsidr)
      ,CHFp_CONST_REAL(r0)
      ,CHFp_CONST_REAL(drr0)
      ,CHFp_CONST_REAL(s_kappa)
      ,CHFp_CONST_REAL(s_delta)
      ,CHFp_FRA1(RBpol)
      ,CHFp_FRA1(dRBpoldt) );

#define FORT_GET_RBPOL_MILLER FORTRAN_NAME( inlineGET_RBPOL_MILLER, inlineGET_RBPOL_MILLER)
#define FORTNT_GET_RBPOL_MILLER FORTRAN_NAME( GET_RBPOL_MILLER, get_rbpol_miller)

inline void FORTRAN_NAME(inlineGET_RBPOL_MILLER, inlineGET_RBPOL_MILLER)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(kappa)
      ,CHFp_CONST_REAL(dpsidr)
      ,CHFp_CONST_REAL(r0)
      ,CHFp_CONST_REAL(drr0)
      ,CHFp_CONST_REAL(s_kappa)
      ,CHFp_CONST_REAL(s_delta)
      ,CHFp_FRA1(RBpol)
      ,CHFp_FRA1(dRBpoldt) )
{
 CH_TIMELEAF("FORT_GET_RBPOL_MILLER");
 FORTRAN_NAME( GET_RBPOL_MILLER ,get_rbpol_miller )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(kappa)
      ,CHFt_CONST_REAL(dpsidr)
      ,CHFt_CONST_REAL(r0)
      ,CHFt_CONST_REAL(drr0)
      ,CHFt_CONST_REAL(s_kappa)
      ,CHFt_CONST_REAL(s_delta)
      ,CHFt_FRA1(RBpol)
      ,CHFt_FRA1(dRBpoldt) );
}
#endif  // GUARDGET_RBPOL_MILLER 

#ifndef GUARDMILLER_DXDXI 
#define GUARDMILLER_DXDXI 
// Prototype for Fortran procedure miller_dxdxi ...
//
void FORTRAN_NAME( MILLER_DXDXI ,miller_dxdxi )(
      CHFp_BOX(box)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(kappa)
      ,CHFp_CONST_FRA(Xi)
      ,CHFp_CONST_INT(destComp)
      ,CHFp_CONST_INT(dirX)
      ,CHFp_CONST_INT(dirXi)
      ,CHFp_FRA(dXdXi) );

#define FORT_MILLER_DXDXI FORTRAN_NAME( inlineMILLER_DXDXI, inlineMILLER_DXDXI)
#define FORTNT_MILLER_DXDXI FORTRAN_NAME( MILLER_DXDXI, miller_dxdxi)

inline void FORTRAN_NAME(inlineMILLER_DXDXI, inlineMILLER_DXDXI)(
      CHFp_BOX(box)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(kappa)
      ,CHFp_CONST_FRA(Xi)
      ,CHFp_CONST_INT(destComp)
      ,CHFp_CONST_INT(dirX)
      ,CHFp_CONST_INT(dirXi)
      ,CHFp_FRA(dXdXi) )
{
 CH_TIMELEAF("FORT_MILLER_DXDXI");
 FORTRAN_NAME( MILLER_DXDXI ,miller_dxdxi )(
      CHFt_BOX(box)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(kappa)
      ,CHFt_CONST_FRA(Xi)
      ,CHFt_CONST_INT(destComp)
      ,CHFt_CONST_INT(dirX)
      ,CHFt_CONST_INT(dirXi)
      ,CHFt_FRA(dXdXi) );
}
#endif  // GUARDMILLER_DXDXI 

#ifndef GUARDGET_MILLER_NC_MAPPED_COORDS 
#define GUARDGET_MILLER_NC_MAPPED_COORDS 
// Prototype for Fortran procedure get_miller_nc_mapped_coords ...
//
void FORTRAN_NAME( GET_MILLER_NC_MAPPED_COORDS ,get_miller_nc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(rmid)
      ,CHFp_CONST_INT(constminorrad)
      ,CHFp_FRA(xi) );

#define FORT_GET_MILLER_NC_MAPPED_COORDS FORTRAN_NAME( inlineGET_MILLER_NC_MAPPED_COORDS, inlineGET_MILLER_NC_MAPPED_COORDS)
#define FORTNT_GET_MILLER_NC_MAPPED_COORDS FORTRAN_NAME( GET_MILLER_NC_MAPPED_COORDS, get_miller_nc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_MILLER_NC_MAPPED_COORDS, inlineGET_MILLER_NC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(rmid)
      ,CHFp_CONST_INT(constminorrad)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_MILLER_NC_MAPPED_COORDS");
 FORTRAN_NAME( GET_MILLER_NC_MAPPED_COORDS ,get_miller_nc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(rmid)
      ,CHFt_CONST_INT(constminorrad)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_MILLER_NC_MAPPED_COORDS 

#ifndef GUARDGET_MILLER_CC_MAPPED_COORDS 
#define GUARDGET_MILLER_CC_MAPPED_COORDS 
// Prototype for Fortran procedure get_miller_cc_mapped_coords ...
//
void FORTRAN_NAME( GET_MILLER_CC_MAPPED_COORDS ,get_miller_cc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(rmid)
      ,CHFp_CONST_INT(constminorrad)
      ,CHFp_FRA(xi) );

#define FORT_GET_MILLER_CC_MAPPED_COORDS FORTRAN_NAME( inlineGET_MILLER_CC_MAPPED_COORDS, inlineGET_MILLER_CC_MAPPED_COORDS)
#define FORTNT_GET_MILLER_CC_MAPPED_COORDS FORTRAN_NAME( GET_MILLER_CC_MAPPED_COORDS, get_miller_cc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_MILLER_CC_MAPPED_COORDS, inlineGET_MILLER_CC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(rmid)
      ,CHFp_CONST_INT(constminorrad)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_MILLER_CC_MAPPED_COORDS");
 FORTRAN_NAME( GET_MILLER_CC_MAPPED_COORDS ,get_miller_cc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(rmid)
      ,CHFt_CONST_INT(constminorrad)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_MILLER_CC_MAPPED_COORDS 

#ifndef GUARDGET_MILLER_FC_MAPPED_COORDS 
#define GUARDGET_MILLER_FC_MAPPED_COORDS 
// Prototype for Fortran procedure get_miller_fc_mapped_coords ...
//
void FORTRAN_NAME( GET_MILLER_FC_MAPPED_COORDS ,get_miller_fc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(rmid)
      ,CHFp_CONST_INT(constminorrad)
      ,CHFp_FRA(xi) );

#define FORT_GET_MILLER_FC_MAPPED_COORDS FORTRAN_NAME( inlineGET_MILLER_FC_MAPPED_COORDS, inlineGET_MILLER_FC_MAPPED_COORDS)
#define FORTNT_GET_MILLER_FC_MAPPED_COORDS FORTRAN_NAME( GET_MILLER_FC_MAPPED_COORDS, get_miller_fc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_MILLER_FC_MAPPED_COORDS, inlineGET_MILLER_FC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(rmid)
      ,CHFp_CONST_INT(constminorrad)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_MILLER_FC_MAPPED_COORDS");
 FORTRAN_NAME( GET_MILLER_FC_MAPPED_COORDS ,get_miller_fc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(rmid)
      ,CHFt_CONST_INT(constminorrad)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_MILLER_FC_MAPPED_COORDS 

}

#endif
