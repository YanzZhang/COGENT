#ifndef _MAXWELLIANKINETICFUNCTIONF_F_H_
#define _MAXWELLIANKINETICFUNCTIONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSET_MAXWELL4D 
#define GUARDSET_MAXWELL4D 
// Prototype for Fortran procedure set_maxwell4d ...
//
void FORTRAN_NAME( SET_MAXWELL4D ,set_maxwell4d )(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coords)
      ,CHFp_CONST_FRA1(dens)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_CONST_FRA1(vshift)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_REAL(mass) );

#define FORT_SET_MAXWELL4D FORTRAN_NAME( inlineSET_MAXWELL4D, inlineSET_MAXWELL4D)
#define FORTNT_SET_MAXWELL4D FORTRAN_NAME( SET_MAXWELL4D, set_maxwell4d)

inline void FORTRAN_NAME(inlineSET_MAXWELL4D, inlineSET_MAXWELL4D)(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coords)
      ,CHFp_CONST_FRA1(dens)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_CONST_FRA1(vshift)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_REAL(mass) )
{
 CH_TIMELEAF("FORT_SET_MAXWELL4D");
 FORTRAN_NAME( SET_MAXWELL4D ,set_maxwell4d )(
      CHFt_FRA(f)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(coords)
      ,CHFt_CONST_FRA1(dens)
      ,CHFt_CONST_FRA1(temp)
      ,CHFt_CONST_FRA1(vshift)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_CONST_REAL(mass) );
}
#endif  // GUARDSET_MAXWELL4D 

#ifndef GUARDSET_CANONICAL_MAXWELL 
#define GUARDSET_CANONICAL_MAXWELL 
// Prototype for Fortran procedure set_canonical_maxwell ...
//
void FORTRAN_NAME( SET_CANONICAL_MAXWELL ,set_canonical_maxwell )(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coords)
      ,CHFp_CONST_FRA(toroidal_coords)
      ,CHFp_CONST_FRA1(B_inj)
      ,CHFp_CONST_FRA1(psi_inj)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_CONST_REAL(psi_p)
      ,CHFp_CONST_REAL(dpsidr_p)
      ,CHFp_CONST_REAL(n_val)
      ,CHFp_CONST_REAL(n_kappa)
      ,CHFp_CONST_REAL(n_width)
      ,CHFp_CONST_REAL(T_val)
      ,CHFp_CONST_REAL(T_kappa)
      ,CHFp_CONST_REAL(T_width)
      ,CHFp_CONST_REALVECT(mode_coeff)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_REAL(larmor) );

#define FORT_SET_CANONICAL_MAXWELL FORTRAN_NAME( inlineSET_CANONICAL_MAXWELL, inlineSET_CANONICAL_MAXWELL)
#define FORTNT_SET_CANONICAL_MAXWELL FORTRAN_NAME( SET_CANONICAL_MAXWELL, set_canonical_maxwell)

inline void FORTRAN_NAME(inlineSET_CANONICAL_MAXWELL, inlineSET_CANONICAL_MAXWELL)(
      CHFp_FRA(f)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(coords)
      ,CHFp_CONST_FRA(toroidal_coords)
      ,CHFp_CONST_FRA1(B_inj)
      ,CHFp_CONST_FRA1(psi_inj)
      ,CHFp_CONST_REAL(RBtor)
      ,CHFp_CONST_REAL(psi_p)
      ,CHFp_CONST_REAL(dpsidr_p)
      ,CHFp_CONST_REAL(n_val)
      ,CHFp_CONST_REAL(n_kappa)
      ,CHFp_CONST_REAL(n_width)
      ,CHFp_CONST_REAL(T_val)
      ,CHFp_CONST_REAL(T_kappa)
      ,CHFp_CONST_REAL(T_width)
      ,CHFp_CONST_REALVECT(mode_coeff)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_CONST_REAL(charge)
      ,CHFp_CONST_REAL(larmor) )
{
 CH_TIMELEAF("FORT_SET_CANONICAL_MAXWELL");
 FORTRAN_NAME( SET_CANONICAL_MAXWELL ,set_canonical_maxwell )(
      CHFt_FRA(f)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(coords)
      ,CHFt_CONST_FRA(toroidal_coords)
      ,CHFt_CONST_FRA1(B_inj)
      ,CHFt_CONST_FRA1(psi_inj)
      ,CHFt_CONST_REAL(RBtor)
      ,CHFt_CONST_REAL(psi_p)
      ,CHFt_CONST_REAL(dpsidr_p)
      ,CHFt_CONST_REAL(n_val)
      ,CHFt_CONST_REAL(n_kappa)
      ,CHFt_CONST_REAL(n_width)
      ,CHFt_CONST_REAL(T_val)
      ,CHFt_CONST_REAL(T_kappa)
      ,CHFt_CONST_REAL(T_width)
      ,CHFt_CONST_REALVECT(mode_coeff)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_CONST_REAL(charge)
      ,CHFt_CONST_REAL(larmor) );
}
#endif  // GUARDSET_CANONICAL_MAXWELL 

}

#endif
