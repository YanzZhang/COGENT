#ifndef _KERNELSF_F_H_
#define _KERNELSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_VEL_CELL 
#define GUARDCOMPUTE_VEL_CELL 
// Prototype for Fortran procedure compute_vel_cell ...
//
void FORTRAN_NAME( COMPUTE_VEL_CELL ,compute_vel_cell )(
      CHFp_FRA(velCell)
      ,CHFp_FRA(velFace)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox) );

#define FORT_COMPUTE_VEL_CELL FORTRAN_NAME( inlineCOMPUTE_VEL_CELL, inlineCOMPUTE_VEL_CELL)
#define FORTNT_COMPUTE_VEL_CELL FORTRAN_NAME( COMPUTE_VEL_CELL, compute_vel_cell)

inline void FORTRAN_NAME(inlineCOMPUTE_VEL_CELL, inlineCOMPUTE_VEL_CELL)(
      CHFp_FRA(velCell)
      ,CHFp_FRA(velFace)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox) )
{
 CH_TIMELEAF("FORT_COMPUTE_VEL_CELL");
 FORTRAN_NAME( COMPUTE_VEL_CELL ,compute_vel_cell )(
      CHFt_FRA(velCell)
      ,CHFt_FRA(velFace)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox) );
}
#endif  // GUARDCOMPUTE_VEL_CELL 

#ifndef GUARDCOMPUTE_PERP_VEL 
#define GUARDCOMPUTE_PERP_VEL 
// Prototype for Fortran procedure compute_perp_vel ...
//
void FORTRAN_NAME( COMPUTE_PERP_VEL ,compute_perp_vel )(
      CHFp_BOX(box)
      ,CHFp_FRA(result)
      ,CHFp_CONST_FRA1(dfn)
      ,CHFp_CONST_FRA(gkVel)
      ,CHFp_CONST_FRA(velCoords)
      ,CHFp_CONST_FRA(b) );

#define FORT_COMPUTE_PERP_VEL FORTRAN_NAME( inlineCOMPUTE_PERP_VEL, inlineCOMPUTE_PERP_VEL)
#define FORTNT_COMPUTE_PERP_VEL FORTRAN_NAME( COMPUTE_PERP_VEL, compute_perp_vel)

inline void FORTRAN_NAME(inlineCOMPUTE_PERP_VEL, inlineCOMPUTE_PERP_VEL)(
      CHFp_BOX(box)
      ,CHFp_FRA(result)
      ,CHFp_CONST_FRA1(dfn)
      ,CHFp_CONST_FRA(gkVel)
      ,CHFp_CONST_FRA(velCoords)
      ,CHFp_CONST_FRA(b) )
{
 CH_TIMELEAF("FORT_COMPUTE_PERP_VEL");
 FORTRAN_NAME( COMPUTE_PERP_VEL ,compute_perp_vel )(
      CHFt_BOX(box)
      ,CHFt_FRA(result)
      ,CHFt_CONST_FRA1(dfn)
      ,CHFt_CONST_FRA(gkVel)
      ,CHFt_CONST_FRA(velCoords)
      ,CHFt_CONST_FRA(b) );
}
#endif  // GUARDCOMPUTE_PERP_VEL 

#ifndef GUARDCOMPUTE_PAR_MOM 
#define GUARDCOMPUTE_PAR_MOM 
// Prototype for Fortran procedure compute_par_mom ...
//
void FORTRAN_NAME( COMPUTE_PAR_MOM ,compute_par_mom )(
      CHFp_FRA(result)
      ,CHFp_CONST_FRA1(parVel)
      ,CHFp_BOX(gridbox) );

#define FORT_COMPUTE_PAR_MOM FORTRAN_NAME( inlineCOMPUTE_PAR_MOM, inlineCOMPUTE_PAR_MOM)
#define FORTNT_COMPUTE_PAR_MOM FORTRAN_NAME( COMPUTE_PAR_MOM, compute_par_mom)

inline void FORTRAN_NAME(inlineCOMPUTE_PAR_MOM, inlineCOMPUTE_PAR_MOM)(
      CHFp_FRA(result)
      ,CHFp_CONST_FRA1(parVel)
      ,CHFp_BOX(gridbox) )
{
 CH_TIMELEAF("FORT_COMPUTE_PAR_MOM");
 FORTRAN_NAME( COMPUTE_PAR_MOM ,compute_par_mom )(
      CHFt_FRA(result)
      ,CHFt_CONST_FRA1(parVel)
      ,CHFt_BOX(gridbox) );
}
#endif  // GUARDCOMPUTE_PAR_MOM 

#ifndef GUARDCOMPUTE_PRESSURE 
#define GUARDCOMPUTE_PRESSURE 
// Prototype for Fortran procedure compute_pressure ...
//
void FORTRAN_NAME( COMPUTE_PRESSURE ,compute_pressure )(
      CHFp_FRA(result)
      ,CHFp_CONST_FRA(velCoords)
      ,CHFp_CONST_FRA1(vparShift)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_BOX(gridbox) );

#define FORT_COMPUTE_PRESSURE FORTRAN_NAME( inlineCOMPUTE_PRESSURE, inlineCOMPUTE_PRESSURE)
#define FORTNT_COMPUTE_PRESSURE FORTRAN_NAME( COMPUTE_PRESSURE, compute_pressure)

inline void FORTRAN_NAME(inlineCOMPUTE_PRESSURE, inlineCOMPUTE_PRESSURE)(
      CHFp_FRA(result)
      ,CHFp_CONST_FRA(velCoords)
      ,CHFp_CONST_FRA1(vparShift)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_BOX(gridbox) )
{
 CH_TIMELEAF("FORT_COMPUTE_PRESSURE");
 FORTRAN_NAME( COMPUTE_PRESSURE ,compute_pressure )(
      CHFt_FRA(result)
      ,CHFt_CONST_FRA(velCoords)
      ,CHFt_CONST_FRA1(vparShift)
      ,CHFt_CONST_FRA1(B)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_BOX(gridbox) );
}
#endif  // GUARDCOMPUTE_PRESSURE 

#ifndef GUARDCOMPUTE_MAXWELLIAN 
#define GUARDCOMPUTE_MAXWELLIAN 
// Prototype for Fortran procedure compute_maxwellian ...
//
void FORTRAN_NAME( COMPUTE_MAXWELLIAN ,compute_maxwellian )(
      CHFp_FRA(result)
      ,CHFp_CONST_FRA(velCoords)
      ,CHFp_CONST_FRA1(vparShift)
      ,CHFp_CONST_FRA1(n)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_BOX(gridbox) );

#define FORT_COMPUTE_MAXWELLIAN FORTRAN_NAME( inlineCOMPUTE_MAXWELLIAN, inlineCOMPUTE_MAXWELLIAN)
#define FORTNT_COMPUTE_MAXWELLIAN FORTRAN_NAME( COMPUTE_MAXWELLIAN, compute_maxwellian)

inline void FORTRAN_NAME(inlineCOMPUTE_MAXWELLIAN, inlineCOMPUTE_MAXWELLIAN)(
      CHFp_FRA(result)
      ,CHFp_CONST_FRA(velCoords)
      ,CHFp_CONST_FRA1(vparShift)
      ,CHFp_CONST_FRA1(n)
      ,CHFp_CONST_FRA1(T)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_CONST_REAL(mass)
      ,CHFp_BOX(gridbox) )
{
 CH_TIMELEAF("FORT_COMPUTE_MAXWELLIAN");
 FORTRAN_NAME( COMPUTE_MAXWELLIAN ,compute_maxwellian )(
      CHFt_FRA(result)
      ,CHFt_CONST_FRA(velCoords)
      ,CHFt_CONST_FRA1(vparShift)
      ,CHFt_CONST_FRA1(n)
      ,CHFt_CONST_FRA1(T)
      ,CHFt_CONST_FRA1(B)
      ,CHFt_CONST_REAL(mass)
      ,CHFt_BOX(gridbox) );
}
#endif  // GUARDCOMPUTE_MAXWELLIAN 

}

#endif
