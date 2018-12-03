#ifndef _NEUTRALSF_F_H_
#define _NEUTRALSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_IONIZATION 
#define GUARDCOMPUTE_IONIZATION 
// Prototype for Fortran procedure compute_ionization ...
//
void FORTRAN_NAME( COMPUTE_IONIZATION ,compute_ionization )(
      CHFp_BOX(box)
      ,CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(fi)
      ,CHFp_CONST_FRA1(nn)
      ,CHFp_CONST_FRA1(sigmaV) );

#define FORT_COMPUTE_IONIZATION FORTRAN_NAME( inlineCOMPUTE_IONIZATION, inlineCOMPUTE_IONIZATION)
#define FORTNT_COMPUTE_IONIZATION FORTRAN_NAME( COMPUTE_IONIZATION, compute_ionization)

inline void FORTRAN_NAME(inlineCOMPUTE_IONIZATION, inlineCOMPUTE_IONIZATION)(
      CHFp_BOX(box)
      ,CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(fi)
      ,CHFp_CONST_FRA1(nn)
      ,CHFp_CONST_FRA1(sigmaV) )
{
 CH_TIMELEAF("FORT_COMPUTE_IONIZATION");
 FORTRAN_NAME( COMPUTE_IONIZATION ,compute_ionization )(
      CHFt_BOX(box)
      ,CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA1(fi)
      ,CHFt_CONST_FRA1(nn)
      ,CHFt_CONST_FRA1(sigmaV) );
}
#endif  // GUARDCOMPUTE_IONIZATION 

#ifndef GUARDCOMPUTE_CHARGE_EXCHANGE 
#define GUARDCOMPUTE_CHARGE_EXCHANGE 
// Prototype for Fortran procedure compute_charge_exchange ...
//
void FORTRAN_NAME( COMPUTE_CHARGE_EXCHANGE ,compute_charge_exchange )(
      CHFp_BOX(box)
      ,CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(fi)
      ,CHFp_CONST_FRA1(ni)
      ,CHFp_CONST_FRA1(Ti)
      ,CHFp_CONST_FRA1(nn)
      ,CHFp_CONST_FRA1(fn) );

#define FORT_COMPUTE_CHARGE_EXCHANGE FORTRAN_NAME( inlineCOMPUTE_CHARGE_EXCHANGE, inlineCOMPUTE_CHARGE_EXCHANGE)
#define FORTNT_COMPUTE_CHARGE_EXCHANGE FORTRAN_NAME( COMPUTE_CHARGE_EXCHANGE, compute_charge_exchange)

inline void FORTRAN_NAME(inlineCOMPUTE_CHARGE_EXCHANGE, inlineCOMPUTE_CHARGE_EXCHANGE)(
      CHFp_BOX(box)
      ,CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(fi)
      ,CHFp_CONST_FRA1(ni)
      ,CHFp_CONST_FRA1(Ti)
      ,CHFp_CONST_FRA1(nn)
      ,CHFp_CONST_FRA1(fn) )
{
 CH_TIMELEAF("FORT_COMPUTE_CHARGE_EXCHANGE");
 FORTRAN_NAME( COMPUTE_CHARGE_EXCHANGE ,compute_charge_exchange )(
      CHFt_BOX(box)
      ,CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA1(fi)
      ,CHFt_CONST_FRA1(ni)
      ,CHFt_CONST_FRA1(Ti)
      ,CHFt_CONST_FRA1(nn)
      ,CHFt_CONST_FRA1(fn) );
}
#endif  // GUARDCOMPUTE_CHARGE_EXCHANGE 

#ifndef GUARDCOMPUTE_MODEL_CHARGE_EXCHANGE 
#define GUARDCOMPUTE_MODEL_CHARGE_EXCHANGE 
// Prototype for Fortran procedure compute_model_charge_exchange ...
//
void FORTRAN_NAME( COMPUTE_MODEL_CHARGE_EXCHANGE ,compute_model_charge_exchange )(
      CHFp_BOX(box)
      ,CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(nn)
      ,CHFp_CONST_FRA1(Ti)
      ,CHFp_CONST_FRA1(u_par)
      ,CHFp_CONST_FRA1(f0)
      ,CHFp_CONST_FRA1(v_par)
      ,CHFp_CONST_REAL(mass) );

#define FORT_COMPUTE_MODEL_CHARGE_EXCHANGE FORTRAN_NAME( inlineCOMPUTE_MODEL_CHARGE_EXCHANGE, inlineCOMPUTE_MODEL_CHARGE_EXCHANGE)
#define FORTNT_COMPUTE_MODEL_CHARGE_EXCHANGE FORTRAN_NAME( COMPUTE_MODEL_CHARGE_EXCHANGE, compute_model_charge_exchange)

inline void FORTRAN_NAME(inlineCOMPUTE_MODEL_CHARGE_EXCHANGE, inlineCOMPUTE_MODEL_CHARGE_EXCHANGE)(
      CHFp_BOX(box)
      ,CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(nn)
      ,CHFp_CONST_FRA1(Ti)
      ,CHFp_CONST_FRA1(u_par)
      ,CHFp_CONST_FRA1(f0)
      ,CHFp_CONST_FRA1(v_par)
      ,CHFp_CONST_REAL(mass) )
{
 CH_TIMELEAF("FORT_COMPUTE_MODEL_CHARGE_EXCHANGE");
 FORTRAN_NAME( COMPUTE_MODEL_CHARGE_EXCHANGE ,compute_model_charge_exchange )(
      CHFt_BOX(box)
      ,CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA1(nn)
      ,CHFt_CONST_FRA1(Ti)
      ,CHFt_CONST_FRA1(u_par)
      ,CHFt_CONST_FRA1(f0)
      ,CHFt_CONST_FRA1(v_par)
      ,CHFt_CONST_REAL(mass) );
}
#endif  // GUARDCOMPUTE_MODEL_CHARGE_EXCHANGE 

}

#endif
