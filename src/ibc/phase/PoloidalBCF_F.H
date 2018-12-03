#ifndef _POLOIDALBCF_F_H_
#define _POLOIDALBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSET_POLOIDAL_BC 
#define GUARDSET_POLOIDAL_BC 
// Prototype for Fortran procedure set_poloidal_bc ...
//
void FORTRAN_NAME( SET_POLOIDAL_BC ,set_poloidal_bc )(
      CHFp_BOX(boundarybox)
      ,CHFp_CONST_FRA1(flipdist)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_R1D(vp,ivphi0)
      ,CHFp_CONST_INT(vp_lo)
      ,CHFp_CONST_R1D(psi,ipsihi0)
      ,CHFp_CONST_INT(psilo)
      ,CHFp_CONST_REAL(psi_lim)
      ,CHFp_FRA1(f)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_INT(ilohi) );

#define FORT_SET_POLOIDAL_BC FORTRAN_NAME( inlineSET_POLOIDAL_BC, inlineSET_POLOIDAL_BC)
#define FORTNT_SET_POLOIDAL_BC FORTRAN_NAME( SET_POLOIDAL_BC, set_poloidal_bc)

inline void FORTRAN_NAME(inlineSET_POLOIDAL_BC, inlineSET_POLOIDAL_BC)(
      CHFp_BOX(boundarybox)
      ,CHFp_CONST_FRA1(flipdist)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_R1D(vp,ivphi0)
      ,CHFp_CONST_INT(vp_lo)
      ,CHFp_CONST_R1D(psi,ipsihi0)
      ,CHFp_CONST_INT(psilo)
      ,CHFp_CONST_REAL(psi_lim)
      ,CHFp_FRA1(f)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_INT(ilohi) )
{
 CH_TIMELEAF("FORT_SET_POLOIDAL_BC");
 FORTRAN_NAME( SET_POLOIDAL_BC ,set_poloidal_bc )(
      CHFt_BOX(boundarybox)
      ,CHFt_CONST_FRA1(flipdist)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_R1D(vp,ivphi0)
      ,CHFt_CONST_INT(vp_lo)
      ,CHFt_CONST_R1D(psi,ipsihi0)
      ,CHFt_CONST_INT(psilo)
      ,CHFt_CONST_REAL(psi_lim)
      ,CHFt_FRA1(f)
      ,CHFt_CONST_REAL(charge)
      ,CHFt_CONST_INT(ilohi) );
}
#endif  // GUARDSET_POLOIDAL_BC 

}

#endif
