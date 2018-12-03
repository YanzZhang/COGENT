#ifndef _MAGGEOMF_F_H_
#define _MAGGEOMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMULT_NJINVERSE 
#define GUARDMULT_NJINVERSE 
// Prototype for Fortran procedure mult_njinverse ...
//
void FORTRAN_NAME( MULT_NJINVERSE ,mult_njinverse )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(efield_mapped)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_FRA(efield) );

#define FORT_MULT_NJINVERSE FORTRAN_NAME( inlineMULT_NJINVERSE, inlineMULT_NJINVERSE)
#define FORTNT_MULT_NJINVERSE FORTRAN_NAME( MULT_NJINVERSE, mult_njinverse)

inline void FORTRAN_NAME(inlineMULT_NJINVERSE, inlineMULT_NJINVERSE)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(efield_mapped)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_FRA(efield) )
{
 CH_TIMELEAF("FORT_MULT_NJINVERSE");
 FORTRAN_NAME( MULT_NJINVERSE ,mult_njinverse )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(efield_mapped)
      ,CHFt_CONST_FRA(njinverse)
      ,CHFt_FRA(efield) );
}
#endif  // GUARDMULT_NJINVERSE 

#ifndef GUARDPROJECT_ONTO_PARALLEL 
#define GUARDPROJECT_ONTO_PARALLEL 
// Prototype for Fortran procedure project_onto_parallel ...
//
void FORTRAN_NAME( PROJECT_ONTO_PARALLEL ,project_onto_parallel )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_FRA(vector) );

#define FORT_PROJECT_ONTO_PARALLEL FORTRAN_NAME( inlinePROJECT_ONTO_PARALLEL, inlinePROJECT_ONTO_PARALLEL)
#define FORTNT_PROJECT_ONTO_PARALLEL FORTRAN_NAME( PROJECT_ONTO_PARALLEL, project_onto_parallel)

inline void FORTRAN_NAME(inlinePROJECT_ONTO_PARALLEL, inlinePROJECT_ONTO_PARALLEL)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_FRA(vector) )
{
 CH_TIMELEAF("FORT_PROJECT_ONTO_PARALLEL");
 FORTRAN_NAME( PROJECT_ONTO_PARALLEL ,project_onto_parallel )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_FRA(vector) );
}
#endif  // GUARDPROJECT_ONTO_PARALLEL 

#ifndef GUARDCOMPUTE_RADIAL_PROJECTION 
#define GUARDCOMPUTE_RADIAL_PROJECTION 
// Prototype for Fortran procedure compute_radial_projection ...
//
void FORTRAN_NAME( COMPUTE_RADIAL_PROJECTION ,compute_radial_projection )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(vector)
      ,CHFp_CONST_FRA1(vector_r) );

#define FORT_COMPUTE_RADIAL_PROJECTION FORTRAN_NAME( inlineCOMPUTE_RADIAL_PROJECTION, inlineCOMPUTE_RADIAL_PROJECTION)
#define FORTNT_COMPUTE_RADIAL_PROJECTION FORTRAN_NAME( COMPUTE_RADIAL_PROJECTION, compute_radial_projection)

inline void FORTRAN_NAME(inlineCOMPUTE_RADIAL_PROJECTION, inlineCOMPUTE_RADIAL_PROJECTION)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(vector)
      ,CHFp_CONST_FRA1(vector_r) )
{
 CH_TIMELEAF("FORT_COMPUTE_RADIAL_PROJECTION");
 FORTRAN_NAME( COMPUTE_RADIAL_PROJECTION ,compute_radial_projection )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_CONST_FRA(vector)
      ,CHFt_CONST_FRA1(vector_r) );
}
#endif  // GUARDCOMPUTE_RADIAL_PROJECTION 

#ifndef GUARDCOMPUTE_POLOIDAL_PROJECTION 
#define GUARDCOMPUTE_POLOIDAL_PROJECTION 
// Prototype for Fortran procedure compute_poloidal_projection ...
//
void FORTRAN_NAME( COMPUTE_POLOIDAL_PROJECTION ,compute_poloidal_projection )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(vector)
      ,CHFp_CONST_FRA1(vector_pol) );

#define FORT_COMPUTE_POLOIDAL_PROJECTION FORTRAN_NAME( inlineCOMPUTE_POLOIDAL_PROJECTION, inlineCOMPUTE_POLOIDAL_PROJECTION)
#define FORTNT_COMPUTE_POLOIDAL_PROJECTION FORTRAN_NAME( COMPUTE_POLOIDAL_PROJECTION, compute_poloidal_projection)

inline void FORTRAN_NAME(inlineCOMPUTE_POLOIDAL_PROJECTION, inlineCOMPUTE_POLOIDAL_PROJECTION)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_CONST_FRA(vector)
      ,CHFp_CONST_FRA1(vector_pol) )
{
 CH_TIMELEAF("FORT_COMPUTE_POLOIDAL_PROJECTION");
 FORTRAN_NAME( COMPUTE_POLOIDAL_PROJECTION ,compute_poloidal_projection )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_CONST_FRA(vector)
      ,CHFt_CONST_FRA1(vector_pol) );
}
#endif  // GUARDCOMPUTE_POLOIDAL_PROJECTION 

#ifndef GUARDCELL_CENTERED_GRAD_COMPONENT 
#define GUARDCELL_CENTERED_GRAD_COMPONENT 
// Prototype for Fortran procedure cell_centered_grad_component ...
//
void FORTRAN_NAME( CELL_CENTERED_GRAD_COMPONENT ,cell_centered_grad_component )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(Efield) );

#define FORT_CELL_CENTERED_GRAD_COMPONENT FORTRAN_NAME( inlineCELL_CENTERED_GRAD_COMPONENT, inlineCELL_CENTERED_GRAD_COMPONENT)
#define FORTNT_CELL_CENTERED_GRAD_COMPONENT FORTRAN_NAME( CELL_CENTERED_GRAD_COMPONENT, cell_centered_grad_component)

inline void FORTRAN_NAME(inlineCELL_CENTERED_GRAD_COMPONENT, inlineCELL_CENTERED_GRAD_COMPONENT)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(Efield) )
{
 CH_TIMELEAF("FORT_CELL_CENTERED_GRAD_COMPONENT");
 FORTRAN_NAME( CELL_CENTERED_GRAD_COMPONENT ,cell_centered_grad_component )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_INT(order)
      ,CHFt_FRA1(Efield) );
}
#endif  // GUARDCELL_CENTERED_GRAD_COMPONENT 

#ifndef GUARDFACE_CENTERED_GRAD_COMPONENT 
#define GUARDFACE_CENTERED_GRAD_COMPONENT 
// Prototype for Fortran procedure face_centered_grad_component ...
//
void FORTRAN_NAME( FACE_CENTERED_GRAD_COMPONENT ,face_centered_grad_component )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(Efield) );

#define FORT_FACE_CENTERED_GRAD_COMPONENT FORTRAN_NAME( inlineFACE_CENTERED_GRAD_COMPONENT, inlineFACE_CENTERED_GRAD_COMPONENT)
#define FORTNT_FACE_CENTERED_GRAD_COMPONENT FORTRAN_NAME( FACE_CENTERED_GRAD_COMPONENT, face_centered_grad_component)

inline void FORTRAN_NAME(inlineFACE_CENTERED_GRAD_COMPONENT, inlineFACE_CENTERED_GRAD_COMPONENT)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(Efield) )
{
 CH_TIMELEAF("FORT_FACE_CENTERED_GRAD_COMPONENT");
 FORTRAN_NAME( FACE_CENTERED_GRAD_COMPONENT ,face_centered_grad_component )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_INT(order)
      ,CHFt_FRA1(Efield) );
}
#endif  // GUARDFACE_CENTERED_GRAD_COMPONENT 

#ifndef GUARDCOMPUTE_FACE_BINVERSE_INTEGRALS 
#define GUARDCOMPUTE_FACE_BINVERSE_INTEGRALS 
// Prototype for Fortran procedure compute_face_Binverse_integrals ...
//
void FORTRAN_NAME( COMPUTE_FACE_BINVERSE_INTEGRALS ,compute_face_binverse_integrals )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_FRA1(integral) );

#define FORT_COMPUTE_FACE_BINVERSE_INTEGRALS FORTRAN_NAME( inlineCOMPUTE_FACE_BINVERSE_INTEGRALS, inlineCOMPUTE_FACE_BINVERSE_INTEGRALS)
#define FORTNT_COMPUTE_FACE_BINVERSE_INTEGRALS FORTRAN_NAME( COMPUTE_FACE_BINVERSE_INTEGRALS, compute_face_binverse_integrals)

inline void FORTRAN_NAME(inlineCOMPUTE_FACE_BINVERSE_INTEGRALS, inlineCOMPUTE_FACE_BINVERSE_INTEGRALS)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_FRA1(integral) )
{
 CH_TIMELEAF("FORT_COMPUTE_FACE_BINVERSE_INTEGRALS");
 FORTRAN_NAME( COMPUTE_FACE_BINVERSE_INTEGRALS ,compute_face_binverse_integrals )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_FRA1(B)
      ,CHFt_FRA1(integral) );
}
#endif  // GUARDCOMPUTE_FACE_BINVERSE_INTEGRALS 

#ifndef GUARDCOMPUTE_FACE_BXGRADB_INTEGRALS 
#define GUARDCOMPUTE_FACE_BXGRADB_INTEGRALS 
// Prototype for Fortran procedure compute_face_bXgradB_integrals ...
//
void FORTRAN_NAME( COMPUTE_FACE_BXGRADB_INTEGRALS ,compute_face_bxgradb_integrals )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_FRA1(integral) );

#define FORT_COMPUTE_FACE_BXGRADB_INTEGRALS FORTRAN_NAME( inlineCOMPUTE_FACE_BXGRADB_INTEGRALS, inlineCOMPUTE_FACE_BXGRADB_INTEGRALS)
#define FORTNT_COMPUTE_FACE_BXGRADB_INTEGRALS FORTRAN_NAME( COMPUTE_FACE_BXGRADB_INTEGRALS, compute_face_bxgradb_integrals)

inline void FORTRAN_NAME(inlineCOMPUTE_FACE_BXGRADB_INTEGRALS, inlineCOMPUTE_FACE_BXGRADB_INTEGRALS)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(B)
      ,CHFp_FRA1(integral) )
{
 CH_TIMELEAF("FORT_COMPUTE_FACE_BXGRADB_INTEGRALS");
 FORTRAN_NAME( COMPUTE_FACE_BXGRADB_INTEGRALS ,compute_face_bxgradb_integrals )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(B)
      ,CHFt_FRA1(integral) );
}
#endif  // GUARDCOMPUTE_FACE_BXGRADB_INTEGRALS 

#ifndef GUARDCOMPUTE_FACE_TRANSVERSE_E 
#define GUARDCOMPUTE_FACE_TRANSVERSE_E 
// Prototype for Fortran procedure compute_face_transverse_E ...
//
void FORTRAN_NAME( COMPUTE_FACE_TRANSVERSE_E ,compute_face_transverse_e )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_FRA1(E) );

#define FORT_COMPUTE_FACE_TRANSVERSE_E FORTRAN_NAME( inlineCOMPUTE_FACE_TRANSVERSE_E, inlineCOMPUTE_FACE_TRANSVERSE_E)
#define FORTNT_COMPUTE_FACE_TRANSVERSE_E FORTRAN_NAME( COMPUTE_FACE_TRANSVERSE_E, compute_face_transverse_e)

inline void FORTRAN_NAME(inlineCOMPUTE_FACE_TRANSVERSE_E, inlineCOMPUTE_FACE_TRANSVERSE_E)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_FRA1(E) )
{
 CH_TIMELEAF("FORT_COMPUTE_FACE_TRANSVERSE_E");
 FORTRAN_NAME( COMPUTE_FACE_TRANSVERSE_E ,compute_face_transverse_e )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_FRA1(E) );
}
#endif  // GUARDCOMPUTE_FACE_TRANSVERSE_E 

#ifndef GUARDCOMPUTE_VOLUME_B_INTEGRALS 
#define GUARDCOMPUTE_VOLUME_B_INTEGRALS 
// Prototype for Fortran procedure compute_volume_B_integrals ...
//
void FORTRAN_NAME( COMPUTE_VOLUME_B_INTEGRALS ,compute_volume_b_integrals )(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(B)
      ,CHFp_CONST_FRA1(jacobian)
      ,CHFp_FRA(integral) );

#define FORT_COMPUTE_VOLUME_B_INTEGRALS FORTRAN_NAME( inlineCOMPUTE_VOLUME_B_INTEGRALS, inlineCOMPUTE_VOLUME_B_INTEGRALS)
#define FORTNT_COMPUTE_VOLUME_B_INTEGRALS FORTRAN_NAME( COMPUTE_VOLUME_B_INTEGRALS, compute_volume_b_integrals)

inline void FORTRAN_NAME(inlineCOMPUTE_VOLUME_B_INTEGRALS, inlineCOMPUTE_VOLUME_B_INTEGRALS)(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(B)
      ,CHFp_CONST_FRA1(jacobian)
      ,CHFp_FRA(integral) )
{
 CH_TIMELEAF("FORT_COMPUTE_VOLUME_B_INTEGRALS");
 FORTRAN_NAME( COMPUTE_VOLUME_B_INTEGRALS ,compute_volume_b_integrals )(
      CHFt_BOX(box)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_FRA(B)
      ,CHFt_CONST_FRA1(jacobian)
      ,CHFt_FRA(integral) );
}
#endif  // GUARDCOMPUTE_VOLUME_B_INTEGRALS 

#ifndef GUARDCOMPUTE_VOLUME_B_DOT_E_INTEGRALS 
#define GUARDCOMPUTE_VOLUME_B_DOT_E_INTEGRALS 
// Prototype for Fortran procedure compute_volume_B_dot_E_integrals ...
//
void FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_E_INTEGRALS ,compute_volume_b_dot_e_integrals )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(E)
      ,CHFp_CONST_FRA(Bintegral)
      ,CHFp_FRA1(integral) );

#define FORT_COMPUTE_VOLUME_B_DOT_E_INTEGRALS FORTRAN_NAME( inlineCOMPUTE_VOLUME_B_DOT_E_INTEGRALS, inlineCOMPUTE_VOLUME_B_DOT_E_INTEGRALS)
#define FORTNT_COMPUTE_VOLUME_B_DOT_E_INTEGRALS FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_E_INTEGRALS, compute_volume_b_dot_e_integrals)

inline void FORTRAN_NAME(inlineCOMPUTE_VOLUME_B_DOT_E_INTEGRALS, inlineCOMPUTE_VOLUME_B_DOT_E_INTEGRALS)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(E)
      ,CHFp_CONST_FRA(Bintegral)
      ,CHFp_FRA1(integral) )
{
 CH_TIMELEAF("FORT_COMPUTE_VOLUME_B_DOT_E_INTEGRALS");
 FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_E_INTEGRALS ,compute_volume_b_dot_e_integrals )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(E)
      ,CHFt_CONST_FRA(Bintegral)
      ,CHFt_FRA1(integral) );
}
#endif  // GUARDCOMPUTE_VOLUME_B_DOT_E_INTEGRALS 

#ifndef GUARDCOMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED 
#define GUARDCOMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED 
// Prototype for Fortran procedure compute_volume_B_dot_E_integrals_field_aligned ...
//
void FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED ,compute_volume_b_dot_e_integrals_field_aligned )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA1(phi_node)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_CONST_FRA(Bintegral)
      ,CHFp_FRA1(integral) );

#define FORT_COMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED FORTRAN_NAME( inlineCOMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED, inlineCOMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED)
#define FORTNT_COMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED, compute_volume_b_dot_e_integrals_field_aligned)

inline void FORTRAN_NAME(inlineCOMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED, inlineCOMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA1(phi_node)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_CONST_FRA(Bintegral)
      ,CHFp_FRA1(integral) )
{
 CH_TIMELEAF("FORT_COMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED");
 FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED ,compute_volume_b_dot_e_integrals_field_aligned )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA1(phi_node)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_FRA(njinverse)
      ,CHFt_CONST_FRA(Bintegral)
      ,CHFt_FRA1(integral) );
}
#endif  // GUARDCOMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED 

#ifndef GUARDCOMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS 
#define GUARDCOMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS 
// Prototype for Fortran procedure compute_volume_B_dot_gradB_integrals ...
//
void FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS ,compute_volume_b_dot_gradb_integrals )(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(B)
      ,CHFp_CONST_FRA(gradB)
      ,CHFp_CONST_FRA1(jacobian)
      ,CHFp_FRA1(integral) );

#define FORT_COMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS FORTRAN_NAME( inlineCOMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS, inlineCOMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS)
#define FORTNT_COMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS, compute_volume_b_dot_gradb_integrals)

inline void FORTRAN_NAME(inlineCOMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS, inlineCOMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS)(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_FRA(B)
      ,CHFp_CONST_FRA(gradB)
      ,CHFp_CONST_FRA1(jacobian)
      ,CHFp_FRA1(integral) )
{
 CH_TIMELEAF("FORT_COMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS");
 FORTRAN_NAME( COMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS ,compute_volume_b_dot_gradb_integrals )(
      CHFt_BOX(box)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_FRA(B)
      ,CHFt_CONST_FRA(gradB)
      ,CHFt_CONST_FRA1(jacobian)
      ,CHFt_FRA1(integral) );
}
#endif  // GUARDCOMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS 

#ifndef GUARDCOMPUTE_ELLIPTIC_OP_COEFF 
#define GUARDCOMPUTE_ELLIPTIC_OP_COEFF 
// Prototype for Fortran procedure compute_elliptic_op_coeff ...
//
void FORTRAN_NAME( COMPUTE_ELLIPTIC_OP_COEFF ,compute_elliptic_op_coeff )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_FRA(perp_coef)
      ,CHFp_FRA(par_coef) );

#define FORT_COMPUTE_ELLIPTIC_OP_COEFF FORTRAN_NAME( inlineCOMPUTE_ELLIPTIC_OP_COEFF, inlineCOMPUTE_ELLIPTIC_OP_COEFF)
#define FORTNT_COMPUTE_ELLIPTIC_OP_COEFF FORTRAN_NAME( COMPUTE_ELLIPTIC_OP_COEFF, compute_elliptic_op_coeff)

inline void FORTRAN_NAME(inlineCOMPUTE_ELLIPTIC_OP_COEFF, inlineCOMPUTE_ELLIPTIC_OP_COEFF)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_FRA(perp_coef)
      ,CHFp_FRA(par_coef) )
{
 CH_TIMELEAF("FORT_COMPUTE_ELLIPTIC_OP_COEFF");
 FORTRAN_NAME( COMPUTE_ELLIPTIC_OP_COEFF ,compute_elliptic_op_coeff )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_FRA(perp_coef)
      ,CHFt_FRA(par_coef) );
}
#endif  // GUARDCOMPUTE_ELLIPTIC_OP_COEFF 

#ifndef GUARDCOMPUTE_ELLIPTIC_OP_COEFF_MAPPED 
#define GUARDCOMPUTE_ELLIPTIC_OP_COEFF_MAPPED 
// Prototype for Fortran procedure compute_elliptic_op_coeff_mapped ...
//
void FORTRAN_NAME( COMPUTE_ELLIPTIC_OP_COEFF_MAPPED ,compute_elliptic_op_coeff_mapped )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(n)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_FRA(coef) );

#define FORT_COMPUTE_ELLIPTIC_OP_COEFF_MAPPED FORTRAN_NAME( inlineCOMPUTE_ELLIPTIC_OP_COEFF_MAPPED, inlineCOMPUTE_ELLIPTIC_OP_COEFF_MAPPED)
#define FORTNT_COMPUTE_ELLIPTIC_OP_COEFF_MAPPED FORTRAN_NAME( COMPUTE_ELLIPTIC_OP_COEFF_MAPPED, compute_elliptic_op_coeff_mapped)

inline void FORTRAN_NAME(inlineCOMPUTE_ELLIPTIC_OP_COEFF_MAPPED, inlineCOMPUTE_ELLIPTIC_OP_COEFF_MAPPED)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(n)
      ,CHFp_CONST_FRA(njinverse)
      ,CHFp_FRA(coef) )
{
 CH_TIMELEAF("FORT_COMPUTE_ELLIPTIC_OP_COEFF_MAPPED");
 FORTRAN_NAME( COMPUTE_ELLIPTIC_OP_COEFF_MAPPED ,compute_elliptic_op_coeff_mapped )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(n)
      ,CHFt_CONST_FRA(njinverse)
      ,CHFt_FRA(coef) );
}
#endif  // GUARDCOMPUTE_ELLIPTIC_OP_COEFF_MAPPED 

#ifndef GUARDCOMPUTE_RADIAL_ELLIPTIC_OP_COEFF 
#define GUARDCOMPUTE_RADIAL_ELLIPTIC_OP_COEFF 
// Prototype for Fortran procedure compute_radial_elliptic_op_coeff ...
//
void FORTRAN_NAME( COMPUTE_RADIAL_ELLIPTIC_OP_COEFF ,compute_radial_elliptic_op_coeff )(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_FRA(rad_coef) );

#define FORT_COMPUTE_RADIAL_ELLIPTIC_OP_COEFF FORTRAN_NAME( inlineCOMPUTE_RADIAL_ELLIPTIC_OP_COEFF, inlineCOMPUTE_RADIAL_ELLIPTIC_OP_COEFF)
#define FORTNT_COMPUTE_RADIAL_ELLIPTIC_OP_COEFF FORTRAN_NAME( COMPUTE_RADIAL_ELLIPTIC_OP_COEFF, compute_radial_elliptic_op_coeff)

inline void FORTRAN_NAME(inlineCOMPUTE_RADIAL_ELLIPTIC_OP_COEFF, inlineCOMPUTE_RADIAL_ELLIPTIC_OP_COEFF)(
      CHFp_BOX(box)
      ,CHFp_CONST_FRA(bunit)
      ,CHFp_FRA(rad_coef) )
{
 CH_TIMELEAF("FORT_COMPUTE_RADIAL_ELLIPTIC_OP_COEFF");
 FORTRAN_NAME( COMPUTE_RADIAL_ELLIPTIC_OP_COEFF ,compute_radial_elliptic_op_coeff )(
      CHFt_BOX(box)
      ,CHFt_CONST_FRA(bunit)
      ,CHFt_FRA(rad_coef) );
}
#endif  // GUARDCOMPUTE_RADIAL_ELLIPTIC_OP_COEFF 

}

#endif
