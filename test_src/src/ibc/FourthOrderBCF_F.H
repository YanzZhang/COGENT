#ifndef _FOURTHORDERBCF_F_H_
#define _FOURTHORDERBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDFOURTH_ORDER_OUTFLOW_BC 
#define GUARDFOURTH_ORDER_OUTFLOW_BC 
// Prototype for Fortran procedure fourth_order_outflow_bc ...
//
void FORTRAN_NAME( FOURTH_ORDER_OUTFLOW_BC ,fourth_order_outflow_bc )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_FRA(vn)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside) );

#define FORT_FOURTH_ORDER_OUTFLOW_BC FORTRAN_NAME( inlineFOURTH_ORDER_OUTFLOW_BC, inlineFOURTH_ORDER_OUTFLOW_BC)
#define FORTNT_FOURTH_ORDER_OUTFLOW_BC FORTRAN_NAME( FOURTH_ORDER_OUTFLOW_BC, fourth_order_outflow_bc)

inline void FORTRAN_NAME(inlineFOURTH_ORDER_OUTFLOW_BC, inlineFOURTH_ORDER_OUTFLOW_BC)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_FRA(vn)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside) )
{
 CH_TIMELEAF("FORT_FOURTH_ORDER_OUTFLOW_BC");
 FORTRAN_NAME( FOURTH_ORDER_OUTFLOW_BC ,fourth_order_outflow_bc )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_FRA(vn)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iside) );
}
#endif  // GUARDFOURTH_ORDER_OUTFLOW_BC 

#ifndef GUARDFOURTH_ORDER_NEUMANN_BC 
#define GUARDFOURTH_ORDER_NEUMANN_BC 
// Prototype for Fortran procedure fourth_order_neumann_bc ...
//
void FORTRAN_NAME( FOURTH_ORDER_NEUMANN_BC ,fourth_order_neumann_bc )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_INTVECT(ng)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside) );

#define FORT_FOURTH_ORDER_NEUMANN_BC FORTRAN_NAME( inlineFOURTH_ORDER_NEUMANN_BC, inlineFOURTH_ORDER_NEUMANN_BC)
#define FORTNT_FOURTH_ORDER_NEUMANN_BC FORTRAN_NAME( FOURTH_ORDER_NEUMANN_BC, fourth_order_neumann_bc)

inline void FORTRAN_NAME(inlineFOURTH_ORDER_NEUMANN_BC, inlineFOURTH_ORDER_NEUMANN_BC)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_INTVECT(ng)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside) )
{
 CH_TIMELEAF("FORT_FOURTH_ORDER_NEUMANN_BC");
 FORTRAN_NAME( FOURTH_ORDER_NEUMANN_BC ,fourth_order_neumann_bc )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_INTVECT(ng)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iside) );
}
#endif  // GUARDFOURTH_ORDER_NEUMANN_BC 

#ifndef GUARDFOURTH_ORDER_DIRICHLET_BC 
#define GUARDFOURTH_ORDER_DIRICHLET_BC 
// Prototype for Fortran procedure fourth_order_dirichlet_bc ...
//
void FORTRAN_NAME( FOURTH_ORDER_DIRICHLET_BC ,fourth_order_dirichlet_bc )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_INTVECT(ng)
      ,CHFp_CONST_FRA(val)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside) );

#define FORT_FOURTH_ORDER_DIRICHLET_BC FORTRAN_NAME( inlineFOURTH_ORDER_DIRICHLET_BC, inlineFOURTH_ORDER_DIRICHLET_BC)
#define FORTNT_FOURTH_ORDER_DIRICHLET_BC FORTRAN_NAME( FOURTH_ORDER_DIRICHLET_BC, fourth_order_dirichlet_bc)

inline void FORTRAN_NAME(inlineFOURTH_ORDER_DIRICHLET_BC, inlineFOURTH_ORDER_DIRICHLET_BC)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_INTVECT(ng)
      ,CHFp_CONST_FRA(val)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside) )
{
 CH_TIMELEAF("FORT_FOURTH_ORDER_DIRICHLET_BC");
 FORTRAN_NAME( FOURTH_ORDER_DIRICHLET_BC ,fourth_order_dirichlet_bc )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_INTVECT(ng)
      ,CHFt_CONST_FRA(val)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iside) );
}
#endif  // GUARDFOURTH_ORDER_DIRICHLET_BC 

#ifndef GUARDSECOND_ORDER_DIRICHLET_BC 
#define GUARDSECOND_ORDER_DIRICHLET_BC 
// Prototype for Fortran procedure second_order_dirichlet_bc ...
//
void FORTRAN_NAME( SECOND_ORDER_DIRICHLET_BC ,second_order_dirichlet_bc )(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside) );

#define FORT_SECOND_ORDER_DIRICHLET_BC FORTRAN_NAME( inlineSECOND_ORDER_DIRICHLET_BC, inlineSECOND_ORDER_DIRICHLET_BC)
#define FORTNT_SECOND_ORDER_DIRICHLET_BC FORTRAN_NAME( SECOND_ORDER_DIRICHLET_BC, second_order_dirichlet_bc)

inline void FORTRAN_NAME(inlineSECOND_ORDER_DIRICHLET_BC, inlineSECOND_ORDER_DIRICHLET_BC)(
      CHFp_FRA(f)
      ,CHFp_BOX(bdrybox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside) )
{
 CH_TIMELEAF("FORT_SECOND_ORDER_DIRICHLET_BC");
 FORTRAN_NAME( SECOND_ORDER_DIRICHLET_BC ,second_order_dirichlet_bc )(
      CHFt_FRA(f)
      ,CHFt_BOX(bdrybox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iside) );
}
#endif  // GUARDSECOND_ORDER_DIRICHLET_BC 

}

#endif
