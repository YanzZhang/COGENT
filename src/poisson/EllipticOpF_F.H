#ifndef _ELLIPTICOPF_F_H_
#define _ELLIPTICOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCELL_CENTERED_FIELD_COMPONENT 
#define GUARDCELL_CENTERED_FIELD_COMPONENT 
// Prototype for Fortran procedure cell_centered_field_component ...
//
void FORTRAN_NAME( CELL_CENTERED_FIELD_COMPONENT ,cell_centered_field_component )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(Efield) );

#define FORT_CELL_CENTERED_FIELD_COMPONENT FORTRAN_NAME( inlineCELL_CENTERED_FIELD_COMPONENT, inlineCELL_CENTERED_FIELD_COMPONENT)
#define FORTNT_CELL_CENTERED_FIELD_COMPONENT FORTRAN_NAME( CELL_CENTERED_FIELD_COMPONENT, cell_centered_field_component)

inline void FORTRAN_NAME(inlineCELL_CENTERED_FIELD_COMPONENT, inlineCELL_CENTERED_FIELD_COMPONENT)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(Efield) )
{
 CH_TIMELEAF("FORT_CELL_CENTERED_FIELD_COMPONENT");
 FORTRAN_NAME( CELL_CENTERED_FIELD_COMPONENT ,cell_centered_field_component )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_INT(order)
      ,CHFt_FRA1(Efield) );
}
#endif  // GUARDCELL_CENTERED_FIELD_COMPONENT 

#ifndef GUARDFACE_CENTERED_FIELD_COMPONENT 
#define GUARDFACE_CENTERED_FIELD_COMPONENT 
// Prototype for Fortran procedure face_centered_field_component ...
//
void FORTRAN_NAME( FACE_CENTERED_FIELD_COMPONENT ,face_centered_field_component )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(Efield) );

#define FORT_FACE_CENTERED_FIELD_COMPONENT FORTRAN_NAME( inlineFACE_CENTERED_FIELD_COMPONENT, inlineFACE_CENTERED_FIELD_COMPONENT)
#define FORTNT_FACE_CENTERED_FIELD_COMPONENT FORTRAN_NAME( FACE_CENTERED_FIELD_COMPONENT, face_centered_field_component)

inline void FORTRAN_NAME(inlineFACE_CENTERED_FIELD_COMPONENT, inlineFACE_CENTERED_FIELD_COMPONENT)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(Efield) )
{
 CH_TIMELEAF("FORT_FACE_CENTERED_FIELD_COMPONENT");
 FORTRAN_NAME( FACE_CENTERED_FIELD_COMPONENT ,face_centered_field_component )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_INT(order)
      ,CHFt_FRA1(Efield) );
}
#endif  // GUARDFACE_CENTERED_FIELD_COMPONENT 

#ifndef GUARDFACE_INTERPOLATE 
#define GUARDFACE_INTERPOLATE 
// Prototype for Fortran procedure face_interpolate ...
//
void FORTRAN_NAME( FACE_INTERPOLATE ,face_interpolate )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(order)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_FRA1(data) );

#define FORT_FACE_INTERPOLATE FORTRAN_NAME( inlineFACE_INTERPOLATE, inlineFACE_INTERPOLATE)
#define FORTNT_FACE_INTERPOLATE FORTRAN_NAME( FACE_INTERPOLATE, face_interpolate)

inline void FORTRAN_NAME(inlineFACE_INTERPOLATE, inlineFACE_INTERPOLATE)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(order)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_FRA1(data) )
{
 CH_TIMELEAF("FORT_FACE_INTERPOLATE");
 FORTRAN_NAME( FACE_INTERPOLATE ,face_interpolate )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(order)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_FRA1(data) );
}
#endif  // GUARDFACE_INTERPOLATE 

#ifndef GUARDEXTRAP_FOR_CC_OPS 
#define GUARDEXTRAP_FOR_CC_OPS 
// Prototype for Fortran procedure extrap_for_cc_ops ...
//
void FORTRAN_NAME( EXTRAP_FOR_CC_OPS ,extrap_for_cc_ops )(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(order)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorbox)
      ,CHFp_FRA(array) );

#define FORT_EXTRAP_FOR_CC_OPS FORTRAN_NAME( inlineEXTRAP_FOR_CC_OPS, inlineEXTRAP_FOR_CC_OPS)
#define FORTNT_EXTRAP_FOR_CC_OPS FORTRAN_NAME( EXTRAP_FOR_CC_OPS, extrap_for_cc_ops)

inline void FORTRAN_NAME(inlineEXTRAP_FOR_CC_OPS, inlineEXTRAP_FOR_CC_OPS)(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(order)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorbox)
      ,CHFp_FRA(array) )
{
 CH_TIMELEAF("FORT_EXTRAP_FOR_CC_OPS");
 FORTRAN_NAME( EXTRAP_FOR_CC_OPS ,extrap_for_cc_ops )(
      CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(order)
      ,CHFt_BOX(facebox)
      ,CHFt_BOX(interiorbox)
      ,CHFt_FRA(array) );
}
#endif  // GUARDEXTRAP_FOR_CC_OPS 

#ifndef GUARDEXTRAP_FOR_FC_OPS 
#define GUARDEXTRAP_FOR_FC_OPS 
// Prototype for Fortran procedure extrap_for_fc_ops ...
//
void FORTRAN_NAME( EXTRAP_FOR_FC_OPS ,extrap_for_fc_ops )(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(order)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorbox)
      ,CHFp_FRA(array) );

#define FORT_EXTRAP_FOR_FC_OPS FORTRAN_NAME( inlineEXTRAP_FOR_FC_OPS, inlineEXTRAP_FOR_FC_OPS)
#define FORTNT_EXTRAP_FOR_FC_OPS FORTRAN_NAME( EXTRAP_FOR_FC_OPS, extrap_for_fc_ops)

inline void FORTRAN_NAME(inlineEXTRAP_FOR_FC_OPS, inlineEXTRAP_FOR_FC_OPS)(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(order)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorbox)
      ,CHFp_FRA(array) )
{
 CH_TIMELEAF("FORT_EXTRAP_FOR_FC_OPS");
 FORTRAN_NAME( EXTRAP_FOR_FC_OPS ,extrap_for_fc_ops )(
      CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(order)
      ,CHFt_BOX(facebox)
      ,CHFt_BOX(interiorbox)
      ,CHFt_FRA(array) );
}
#endif  // GUARDEXTRAP_FOR_FC_OPS 

}

#endif
