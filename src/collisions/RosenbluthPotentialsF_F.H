#ifndef _ROSENBLUTHPOTENTIALSF_F_H_
#define _ROSENBLUTHPOTENTIALSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDIMPOSE_MULTIPOLE_BC1 
#define GUARDIMPOSE_MULTIPOLE_BC1 
// Prototype for Fortran procedure impose_multipole_bc1 ...
//
void FORTRAN_NAME( IMPOSE_MULTIPOLE_BC1 ,impose_multipole_bc1 )(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA(mcoeff)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(p) );

#define FORT_IMPOSE_MULTIPOLE_BC1 FORTRAN_NAME( inlineIMPOSE_MULTIPOLE_BC1, inlineIMPOSE_MULTIPOLE_BC1)
#define FORTNT_IMPOSE_MULTIPOLE_BC1 FORTRAN_NAME( IMPOSE_MULTIPOLE_BC1, impose_multipole_bc1)

inline void FORTRAN_NAME(inlineIMPOSE_MULTIPOLE_BC1, inlineIMPOSE_MULTIPOLE_BC1)(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA(mcoeff)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_IMPOSE_MULTIPOLE_BC1");
 FORTRAN_NAME( IMPOSE_MULTIPOLE_BC1 ,impose_multipole_bc1 )(
      CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA(mcoeff)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDIMPOSE_MULTIPOLE_BC1 

#ifndef GUARDIMPOSE_MULTIPOLE_BC2 
#define GUARDIMPOSE_MULTIPOLE_BC2 
// Prototype for Fortran procedure impose_multipole_bc2 ...
//
void FORTRAN_NAME( IMPOSE_MULTIPOLE_BC2 ,impose_multipole_bc2 )(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA(mcoeff_rho)
      ,CHFp_CONST_FRA(mcoeff_phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(p) );

#define FORT_IMPOSE_MULTIPOLE_BC2 FORTRAN_NAME( inlineIMPOSE_MULTIPOLE_BC2, inlineIMPOSE_MULTIPOLE_BC2)
#define FORTNT_IMPOSE_MULTIPOLE_BC2 FORTRAN_NAME( IMPOSE_MULTIPOLE_BC2, impose_multipole_bc2)

inline void FORTRAN_NAME(inlineIMPOSE_MULTIPOLE_BC2, inlineIMPOSE_MULTIPOLE_BC2)(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA(mcoeff_rho)
      ,CHFp_CONST_FRA(mcoeff_phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_IMPOSE_MULTIPOLE_BC2");
 FORTRAN_NAME( IMPOSE_MULTIPOLE_BC2 ,impose_multipole_bc2 )(
      CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA(mcoeff_rho)
      ,CHFt_CONST_FRA(mcoeff_phi)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDIMPOSE_MULTIPOLE_BC2 

#ifndef GUARDSEPARATE_PHI1_INNER 
#define GUARDSEPARATE_PHI1_INNER 
// Prototype for Fortran procedure separate_phi1_inner ...
//
void FORTRAN_NAME( SEPARATE_PHI1_INNER ,separate_phi1_inner )(
      CHFp_FRA1(phi_inner)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(p) );

#define FORT_SEPARATE_PHI1_INNER FORTRAN_NAME( inlineSEPARATE_PHI1_INNER, inlineSEPARATE_PHI1_INNER)
#define FORTNT_SEPARATE_PHI1_INNER FORTRAN_NAME( SEPARATE_PHI1_INNER, separate_phi1_inner)

inline void FORTRAN_NAME(inlineSEPARATE_PHI1_INNER, inlineSEPARATE_PHI1_INNER)(
      CHFp_FRA1(phi_inner)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(Nvpar)
      ,CHFp_CONST_INT(Nmu)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_SEPARATE_PHI1_INNER");
 FORTRAN_NAME( SEPARATE_PHI1_INNER ,separate_phi1_inner )(
      CHFt_FRA1(phi_inner)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(Nvpar)
      ,CHFt_CONST_INT(Nmu)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDSEPARATE_PHI1_INNER 

#ifndef GUARDCOMPUTE_MULT_KERNELS 
#define GUARDCOMPUTE_MULT_KERNELS 
// Prototype for Fortran procedure compute_mult_kernels ...
//
void FORTRAN_NAME( COMPUTE_MULT_KERNELS ,compute_mult_kernels )(
      CHFp_FRA(kernel)
      ,CHFp_CONST_FRA1(rho)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(p) );

#define FORT_COMPUTE_MULT_KERNELS FORTRAN_NAME( inlineCOMPUTE_MULT_KERNELS, inlineCOMPUTE_MULT_KERNELS)
#define FORTNT_COMPUTE_MULT_KERNELS FORTRAN_NAME( COMPUTE_MULT_KERNELS, compute_mult_kernels)

inline void FORTRAN_NAME(inlineCOMPUTE_MULT_KERNELS, inlineCOMPUTE_MULT_KERNELS)(
      CHFp_FRA(kernel)
      ,CHFp_CONST_FRA1(rho)
      ,CHFp_CONST_FRA1(bmag)
      ,CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_COMPUTE_MULT_KERNELS");
 FORTRAN_NAME( COMPUTE_MULT_KERNELS ,compute_mult_kernels )(
      CHFt_FRA(kernel)
      ,CHFt_CONST_FRA1(rho)
      ,CHFt_CONST_FRA1(bmag)
      ,CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDCOMPUTE_MULT_KERNELS 

}

#endif
