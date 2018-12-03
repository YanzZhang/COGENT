#ifndef _SLABBLOCKCOORDSYSF_F_H_
#define _SLABBLOCKCOORDSYSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGET_SLAB_FIELD_DATA 
#define GUARDGET_SLAB_FIELD_DATA 
// Prototype for Fortran procedure get_slab_field_data ...
//
void FORTRAN_NAME( GET_SLAB_FIELD_DATA ,get_slab_field_data )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(xix)
      ,CHFp_CONST_REAL(ByInner)
      ,CHFp_CONST_REAL(ByOuter)
      ,CHFp_CONST_REAL(BzInner)
      ,CHFp_CONST_REAL(BzOuter)
      ,CHFp_CONST_REAL(xmax)
      ,CHFp_CONST_REAL(ximax)
      ,CHFp_FRA(b_pt)
      ,CHFp_FRA1(Bmag_pt)
      ,CHFp_FRA(bunit_pt)
      ,CHFp_FRA(gradb_pt)
      ,CHFp_FRA(curlb_pt)
      ,CHFp_FRA1(bdotcurlb_pt) );

#define FORT_GET_SLAB_FIELD_DATA FORTRAN_NAME( inlineGET_SLAB_FIELD_DATA, inlineGET_SLAB_FIELD_DATA)
#define FORTNT_GET_SLAB_FIELD_DATA FORTRAN_NAME( GET_SLAB_FIELD_DATA, get_slab_field_data)

inline void FORTRAN_NAME(inlineGET_SLAB_FIELD_DATA, inlineGET_SLAB_FIELD_DATA)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(xix)
      ,CHFp_CONST_REAL(ByInner)
      ,CHFp_CONST_REAL(ByOuter)
      ,CHFp_CONST_REAL(BzInner)
      ,CHFp_CONST_REAL(BzOuter)
      ,CHFp_CONST_REAL(xmax)
      ,CHFp_CONST_REAL(ximax)
      ,CHFp_FRA(b_pt)
      ,CHFp_FRA1(Bmag_pt)
      ,CHFp_FRA(bunit_pt)
      ,CHFp_FRA(gradb_pt)
      ,CHFp_FRA(curlb_pt)
      ,CHFp_FRA1(bdotcurlb_pt) )
{
 CH_TIMELEAF("FORT_GET_SLAB_FIELD_DATA");
 FORTRAN_NAME( GET_SLAB_FIELD_DATA ,get_slab_field_data )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(xix)
      ,CHFt_CONST_REAL(ByInner)
      ,CHFt_CONST_REAL(ByOuter)
      ,CHFt_CONST_REAL(BzInner)
      ,CHFt_CONST_REAL(BzOuter)
      ,CHFt_CONST_REAL(xmax)
      ,CHFt_CONST_REAL(ximax)
      ,CHFt_FRA(b_pt)
      ,CHFt_FRA1(Bmag_pt)
      ,CHFt_FRA(bunit_pt)
      ,CHFt_FRA(gradb_pt)
      ,CHFt_FRA(curlb_pt)
      ,CHFt_FRA1(bdotcurlb_pt) );
}
#endif  // GUARDGET_SLAB_FIELD_DATA 

}

#endif
