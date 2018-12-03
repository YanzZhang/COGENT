#ifndef _SINGLENULLPHASEBCF_F_H_
#define _SINGLENULLPHASEBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSET_LOGICAL_SHEATH_BC 
#define GUARDSET_LOGICAL_SHEATH_BC 
// Prototype for Fortran procedure set_logical_sheath_bc ...
//
void FORTRAN_NAME( SET_LOGICAL_SHEATH_BC ,set_logical_sheath_bc )(
      CHFp_FRA1(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_FRA1(f_rflct)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA1(vn)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_INT(iside) );

#define FORT_SET_LOGICAL_SHEATH_BC FORTRAN_NAME( inlineSET_LOGICAL_SHEATH_BC, inlineSET_LOGICAL_SHEATH_BC)
#define FORTNT_SET_LOGICAL_SHEATH_BC FORTRAN_NAME( SET_LOGICAL_SHEATH_BC, set_logical_sheath_bc)

inline void FORTRAN_NAME(inlineSET_LOGICAL_SHEATH_BC, inlineSET_LOGICAL_SHEATH_BC)(
      CHFp_FRA1(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_FRA1(f_rflct)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA1(vn)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_INT(iside) )
{
 CH_TIMELEAF("FORT_SET_LOGICAL_SHEATH_BC");
 FORTRAN_NAME( SET_LOGICAL_SHEATH_BC ,set_logical_sheath_bc )(
      CHFt_FRA1(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_FRA1(f_rflct)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_FRA1(vn)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_REAL(charge)
      ,CHFt_CONST_INT(iside) );
}
#endif  // GUARDSET_LOGICAL_SHEATH_BC 

#ifndef GUARDSET_POLOIDAL_DIVERTER_BC 
#define GUARDSET_POLOIDAL_DIVERTER_BC 
// Prototype for Fortran procedure set_poloidal_diverter_bc ...
//
void FORTRAN_NAME( SET_POLOIDAL_DIVERTER_BC ,set_poloidal_diverter_bc )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_FRA(f_rflct)
      ,CHFp_CONST_FRA(vn)
      ,CHFp_CONST_FRA1(qphi)
      ,CHFp_CONST_INT(iside) );

#define FORT_SET_POLOIDAL_DIVERTER_BC FORTRAN_NAME( inlineSET_POLOIDAL_DIVERTER_BC, inlineSET_POLOIDAL_DIVERTER_BC)
#define FORTNT_SET_POLOIDAL_DIVERTER_BC FORTRAN_NAME( SET_POLOIDAL_DIVERTER_BC, set_poloidal_diverter_bc)

inline void FORTRAN_NAME(inlineSET_POLOIDAL_DIVERTER_BC, inlineSET_POLOIDAL_DIVERTER_BC)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_FRA(f_rflct)
      ,CHFp_CONST_FRA(vn)
      ,CHFp_CONST_FRA1(qphi)
      ,CHFp_CONST_INT(iside) )
{
 CH_TIMELEAF("FORT_SET_POLOIDAL_DIVERTER_BC");
 FORTRAN_NAME( SET_POLOIDAL_DIVERTER_BC ,set_poloidal_diverter_bc )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_FRA(f_rflct)
      ,CHFt_CONST_FRA(vn)
      ,CHFt_CONST_FRA1(qphi)
      ,CHFt_CONST_INT(iside) );
}
#endif  // GUARDSET_POLOIDAL_DIVERTER_BC 

}

#endif
