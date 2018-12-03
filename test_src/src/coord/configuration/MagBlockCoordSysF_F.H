#ifndef _MAGBLOCKCOORDSYSF_F_H_
#define _MAGBLOCKCOORDSYSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGET_NC_MAPPED_COORDS 
#define GUARDGET_NC_MAPPED_COORDS 
// Prototype for Fortran procedure get_nc_mapped_coords ...
//
void FORTRAN_NAME( GET_NC_MAPPED_COORDS ,get_nc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) );

#define FORT_GET_NC_MAPPED_COORDS FORTRAN_NAME( inlineGET_NC_MAPPED_COORDS, inlineGET_NC_MAPPED_COORDS)
#define FORTNT_GET_NC_MAPPED_COORDS FORTRAN_NAME( GET_NC_MAPPED_COORDS, get_nc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_NC_MAPPED_COORDS, inlineGET_NC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_NC_MAPPED_COORDS");
 FORTRAN_NAME( GET_NC_MAPPED_COORDS ,get_nc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_NC_MAPPED_COORDS 

#ifndef GUARDGET_CC_MAPPED_COORDS 
#define GUARDGET_CC_MAPPED_COORDS 
// Prototype for Fortran procedure get_cc_mapped_coords ...
//
void FORTRAN_NAME( GET_CC_MAPPED_COORDS ,get_cc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) );

#define FORT_GET_CC_MAPPED_COORDS FORTRAN_NAME( inlineGET_CC_MAPPED_COORDS, inlineGET_CC_MAPPED_COORDS)
#define FORTNT_GET_CC_MAPPED_COORDS FORTRAN_NAME( GET_CC_MAPPED_COORDS, get_cc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_CC_MAPPED_COORDS, inlineGET_CC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_CC_MAPPED_COORDS");
 FORTRAN_NAME( GET_CC_MAPPED_COORDS ,get_cc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_CC_MAPPED_COORDS 

#ifndef GUARDGET_FC_MAPPED_COORDS 
#define GUARDGET_FC_MAPPED_COORDS 
// Prototype for Fortran procedure get_fc_mapped_coords ...
//
void FORTRAN_NAME( GET_FC_MAPPED_COORDS ,get_fc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) );

#define FORT_GET_FC_MAPPED_COORDS FORTRAN_NAME( inlineGET_FC_MAPPED_COORDS, inlineGET_FC_MAPPED_COORDS)
#define FORTNT_GET_FC_MAPPED_COORDS FORTRAN_NAME( GET_FC_MAPPED_COORDS, get_fc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_FC_MAPPED_COORDS, inlineGET_FC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_FC_MAPPED_COORDS");
 FORTRAN_NAME( GET_FC_MAPPED_COORDS ,get_fc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_FC_MAPPED_COORDS 

#ifndef GUARDINCREMENTLAPLACIAN2 
#define GUARDINCREMENTLAPLACIAN2 
// Prototype for Fortran procedure INCREMENTLAPLACIAN2 ...
//
void FORTRAN_NAME( INCREMENTLAPLACIAN2 ,incrementlaplacian2 )(
      CHFp_FRA(lapPhi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(factor) );

#define FORT_INCREMENTLAPLACIAN2 FORTRAN_NAME( inlineINCREMENTLAPLACIAN2, inlineINCREMENTLAPLACIAN2)
#define FORTNT_INCREMENTLAPLACIAN2 FORTRAN_NAME( INCREMENTLAPLACIAN2, incrementlaplacian2)

inline void FORTRAN_NAME(inlineINCREMENTLAPLACIAN2, inlineINCREMENTLAPLACIAN2)(
      CHFp_FRA(lapPhi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(factor) )
{
 CH_TIMELEAF("FORT_INCREMENTLAPLACIAN2");
 FORTRAN_NAME( INCREMENTLAPLACIAN2 ,incrementlaplacian2 )(
      CHFt_FRA(lapPhi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(gridBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(factor) );
}
#endif  // GUARDINCREMENTLAPLACIAN2 

#ifndef GUARDUNIT_FS_TANGENT 
#define GUARDUNIT_FS_TANGENT 
// Prototype for Fortran procedure unit_fs_tangent ...
//
void FORTRAN_NAME( UNIT_FS_TANGENT ,unit_fs_tangent )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(dXdxi)
      ,CHFp_FRA(data) );

#define FORT_UNIT_FS_TANGENT FORTRAN_NAME( inlineUNIT_FS_TANGENT, inlineUNIT_FS_TANGENT)
#define FORTNT_UNIT_FS_TANGENT FORTRAN_NAME( UNIT_FS_TANGENT, unit_fs_tangent)

inline void FORTRAN_NAME(inlineUNIT_FS_TANGENT, inlineUNIT_FS_TANGENT)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(dXdxi)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_UNIT_FS_TANGENT");
 FORTRAN_NAME( UNIT_FS_TANGENT ,unit_fs_tangent )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(dXdxi)
      ,CHFt_FRA(data) );
}
#endif  // GUARDUNIT_FS_TANGENT 

#ifndef GUARDUNIT_FS_NORMAL 
#define GUARDUNIT_FS_NORMAL 
// Prototype for Fortran procedure unit_fs_normal ...
//
void FORTRAN_NAME( UNIT_FS_NORMAL ,unit_fs_normal )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(dXdxi)
      ,CHFp_FRA(data) );

#define FORT_UNIT_FS_NORMAL FORTRAN_NAME( inlineUNIT_FS_NORMAL, inlineUNIT_FS_NORMAL)
#define FORTNT_UNIT_FS_NORMAL FORTRAN_NAME( UNIT_FS_NORMAL, unit_fs_normal)

inline void FORTRAN_NAME(inlineUNIT_FS_NORMAL, inlineUNIT_FS_NORMAL)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(dXdxi)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_UNIT_FS_NORMAL");
 FORTRAN_NAME( UNIT_FS_NORMAL ,unit_fs_normal )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(dXdxi)
      ,CHFt_FRA(data) );
}
#endif  // GUARDUNIT_FS_NORMAL 

#ifndef GUARDUNIT_RADIAL 
#define GUARDUNIT_RADIAL 
// Prototype for Fortran procedure unit_radial ...
//
void FORTRAN_NAME( UNIT_RADIAL ,unit_radial )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(dXdxi)
      ,CHFp_FRA(data) );

#define FORT_UNIT_RADIAL FORTRAN_NAME( inlineUNIT_RADIAL, inlineUNIT_RADIAL)
#define FORTNT_UNIT_RADIAL FORTRAN_NAME( UNIT_RADIAL, unit_radial)

inline void FORTRAN_NAME(inlineUNIT_RADIAL, inlineUNIT_RADIAL)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(dXdxi)
      ,CHFp_FRA(data) )
{
 CH_TIMELEAF("FORT_UNIT_RADIAL");
 FORTRAN_NAME( UNIT_RADIAL ,unit_radial )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(dXdxi)
      ,CHFt_FRA(data) );
}
#endif  // GUARDUNIT_RADIAL 

#ifndef GUARDGRADF_FACTOR 
#define GUARDGRADF_FACTOR 
// Prototype for Fortran procedure gradf_factor ...
//
void FORTRAN_NAME( GRADF_FACTOR ,gradf_factor )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(dXdxi)
      ,CHFp_FRA1(data) );

#define FORT_GRADF_FACTOR FORTRAN_NAME( inlineGRADF_FACTOR, inlineGRADF_FACTOR)
#define FORTNT_GRADF_FACTOR FORTRAN_NAME( GRADF_FACTOR, gradf_factor)

inline void FORTRAN_NAME(inlineGRADF_FACTOR, inlineGRADF_FACTOR)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(dXdxi)
      ,CHFp_FRA1(data) )
{
 CH_TIMELEAF("FORT_GRADF_FACTOR");
 FORTRAN_NAME( GRADF_FACTOR ,gradf_factor )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(dXdxi)
      ,CHFt_FRA1(data) );
}
#endif  // GUARDGRADF_FACTOR 

#ifndef GUARDMAG_BLOCK_PROJECT 
#define GUARDMAG_BLOCK_PROJECT 
// Prototype for Fortran procedure mag_block_project ...
//
void FORTRAN_NAME( MAG_BLOCK_PROJECT ,mag_block_project )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(vec_src)
      ,CHFp_FRA(vec_dst) );

#define FORT_MAG_BLOCK_PROJECT FORTRAN_NAME( inlineMAG_BLOCK_PROJECT, inlineMAG_BLOCK_PROJECT)
#define FORTNT_MAG_BLOCK_PROJECT FORTRAN_NAME( MAG_BLOCK_PROJECT, mag_block_project)

inline void FORTRAN_NAME(inlineMAG_BLOCK_PROJECT, inlineMAG_BLOCK_PROJECT)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(vec_src)
      ,CHFp_FRA(vec_dst) )
{
 CH_TIMELEAF("FORT_MAG_BLOCK_PROJECT");
 FORTRAN_NAME( MAG_BLOCK_PROJECT ,mag_block_project )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(vec_src)
      ,CHFt_FRA(vec_dst) );
}
#endif  // GUARDMAG_BLOCK_PROJECT 

#ifndef GUARDMAG_BLOCK_PSITHETA_PROJECTIONS 
#define GUARDMAG_BLOCK_PSITHETA_PROJECTIONS 
// Prototype for Fortran procedure mag_block_psitheta_projections ...
//
void FORTRAN_NAME( MAG_BLOCK_PSITHETA_PROJECTIONS ,mag_block_psitheta_projections )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(vec_psi)
      ,CHFp_CONST_FRA(vec_theta)
      ,CHFp_FRA(vec_dst) );

#define FORT_MAG_BLOCK_PSITHETA_PROJECTIONS FORTRAN_NAME( inlineMAG_BLOCK_PSITHETA_PROJECTIONS, inlineMAG_BLOCK_PSITHETA_PROJECTIONS)
#define FORTNT_MAG_BLOCK_PSITHETA_PROJECTIONS FORTRAN_NAME( MAG_BLOCK_PSITHETA_PROJECTIONS, mag_block_psitheta_projections)

inline void FORTRAN_NAME(inlineMAG_BLOCK_PSITHETA_PROJECTIONS, inlineMAG_BLOCK_PSITHETA_PROJECTIONS)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(vec_psi)
      ,CHFp_CONST_FRA(vec_theta)
      ,CHFp_FRA(vec_dst) )
{
 CH_TIMELEAF("FORT_MAG_BLOCK_PSITHETA_PROJECTIONS");
 FORTRAN_NAME( MAG_BLOCK_PSITHETA_PROJECTIONS ,mag_block_psitheta_projections )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(vec_psi)
      ,CHFt_CONST_FRA(vec_theta)
      ,CHFt_FRA(vec_dst) );
}
#endif  // GUARDMAG_BLOCK_PSITHETA_PROJECTIONS 

#ifndef GUARDGET_NODAL_FIELD_DATA 
#define GUARDGET_NODAL_FIELD_DATA 
// Prototype for Fortran procedure get_nodal_field_data ...
//
void FORTRAN_NAME( GET_NODAL_FIELD_DATA ,get_nodal_field_data )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(RZ)
      ,CHFp_CONST_FRA1(psi)
      ,CHFp_CONST_FRA(RB)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_FRA(A)
      ,CHFp_FRA(bunit)
      ,CHFp_FRA(Bmag) );

#define FORT_GET_NODAL_FIELD_DATA FORTRAN_NAME( inlineGET_NODAL_FIELD_DATA, inlineGET_NODAL_FIELD_DATA)
#define FORTNT_GET_NODAL_FIELD_DATA FORTRAN_NAME( GET_NODAL_FIELD_DATA, get_nodal_field_data)

inline void FORTRAN_NAME(inlineGET_NODAL_FIELD_DATA, inlineGET_NODAL_FIELD_DATA)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(RZ)
      ,CHFp_CONST_FRA1(psi)
      ,CHFp_CONST_FRA(RB)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_FRA(A)
      ,CHFp_FRA(bunit)
      ,CHFp_FRA(Bmag) )
{
 CH_TIMELEAF("FORT_GET_NODAL_FIELD_DATA");
 FORTRAN_NAME( GET_NODAL_FIELD_DATA ,get_nodal_field_data )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(RZ)
      ,CHFt_CONST_FRA1(psi)
      ,CHFt_CONST_FRA(RB)
      ,CHFt_CONST_REAL(RBtor)
      ,CHFt_FRA(A)
      ,CHFt_FRA(bunit)
      ,CHFt_FRA(Bmag) );
}
#endif  // GUARDGET_NODAL_FIELD_DATA 

#ifndef GUARDGET_FIELD_DATA 
#define GUARDGET_FIELD_DATA 
// Prototype for Fortran procedure get_field_data ...
//
void FORTRAN_NAME( GET_FIELD_DATA ,get_field_data )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(RZ)
      ,CHFp_CONST_FRA(RB)
      ,CHFp_CONST_FRA(dRBdR)
      ,CHFp_CONST_FRA(dRBdZ)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_FRA(B)
      ,CHFp_FRA1(Bmag)
      ,CHFp_FRA(bunit)
      ,CHFp_FRA(gradB)
      ,CHFp_FRA(curlb)
      ,CHFp_FRA1(bdotcurlb) );

#define FORT_GET_FIELD_DATA FORTRAN_NAME( inlineGET_FIELD_DATA, inlineGET_FIELD_DATA)
#define FORTNT_GET_FIELD_DATA FORTRAN_NAME( GET_FIELD_DATA, get_field_data)

inline void FORTRAN_NAME(inlineGET_FIELD_DATA, inlineGET_FIELD_DATA)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(RZ)
      ,CHFp_CONST_FRA(RB)
      ,CHFp_CONST_FRA(dRBdR)
      ,CHFp_CONST_FRA(dRBdZ)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_FRA(B)
      ,CHFp_FRA1(Bmag)
      ,CHFp_FRA(bunit)
      ,CHFp_FRA(gradB)
      ,CHFp_FRA(curlb)
      ,CHFp_FRA1(bdotcurlb) )
{
 CH_TIMELEAF("FORT_GET_FIELD_DATA");
 FORTRAN_NAME( GET_FIELD_DATA ,get_field_data )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(RZ)
      ,CHFt_CONST_FRA(RB)
      ,CHFt_CONST_FRA(dRBdR)
      ,CHFt_CONST_FRA(dRBdZ)
      ,CHFt_CONST_REAL(RBtor)
      ,CHFt_FRA(B)
      ,CHFt_FRA1(Bmag)
      ,CHFt_FRA(bunit)
      ,CHFt_FRA(gradB)
      ,CHFt_FRA(curlb)
      ,CHFt_FRA1(bdotcurlb) );
}
#endif  // GUARDGET_FIELD_DATA 

}

#endif
