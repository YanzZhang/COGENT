#ifndef _POSITIVITYPOSTPROCESSORF_F_H_
#define _POSITIVITYPOSTPROCESSORF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDREDISTRIBUTENEGATIVES 
#define GUARDREDISTRIBUTENEGATIVES 
// Prototype for Fortran procedure REDISTRIBUTENEGATIVES ...
//
void FORTRAN_NAME( REDISTRIBUTENEGATIVES ,redistributenegatives )(
      CHFp_FRA(phi)
      ,CHFp_BOX(interiorBox)
      ,CHFp_BOX(neighborBox)
      ,CHFp_CONST_REAL(refVal)
      ,CHFp_INT(unable) );

#define FORT_REDISTRIBUTENEGATIVES FORTRAN_NAME( inlineREDISTRIBUTENEGATIVES, inlineREDISTRIBUTENEGATIVES)
#define FORTNT_REDISTRIBUTENEGATIVES FORTRAN_NAME( REDISTRIBUTENEGATIVES, redistributenegatives)

inline void FORTRAN_NAME(inlineREDISTRIBUTENEGATIVES, inlineREDISTRIBUTENEGATIVES)(
      CHFp_FRA(phi)
      ,CHFp_BOX(interiorBox)
      ,CHFp_BOX(neighborBox)
      ,CHFp_CONST_REAL(refVal)
      ,CHFp_INT(unable) )
{
 CH_TIMELEAF("FORT_REDISTRIBUTENEGATIVES");
 FORTRAN_NAME( REDISTRIBUTENEGATIVES ,redistributenegatives )(
      CHFt_FRA(phi)
      ,CHFt_BOX(interiorBox)
      ,CHFt_BOX(neighborBox)
      ,CHFt_CONST_REAL(refVal)
      ,CHFt_INT(unable) );
}
#endif  // GUARDREDISTRIBUTENEGATIVES 

#ifndef GUARDFINDANYNEGATIVES 
#define GUARDFINDANYNEGATIVES 
// Prototype for Fortran procedure FINDANYNEGATIVES ...
//
void FORTRAN_NAME( FINDANYNEGATIVES ,findanynegatives )(
      CHFp_INT(count)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(interiorBox) );

#define FORT_FINDANYNEGATIVES FORTRAN_NAME( inlineFINDANYNEGATIVES, inlineFINDANYNEGATIVES)
#define FORTNT_FINDANYNEGATIVES FORTRAN_NAME( FINDANYNEGATIVES, findanynegatives)

inline void FORTRAN_NAME(inlineFINDANYNEGATIVES, inlineFINDANYNEGATIVES)(
      CHFp_INT(count)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(interiorBox) )
{
 CH_TIMELEAF("FORT_FINDANYNEGATIVES");
 FORTRAN_NAME( FINDANYNEGATIVES ,findanynegatives )(
      CHFt_INT(count)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(interiorBox) );
}
#endif  // GUARDFINDANYNEGATIVES 

#ifndef GUARDFINDALLNEGATIVES 
#define GUARDFINDALLNEGATIVES 
// Prototype for Fortran procedure FINDALLNEGATIVES ...
//
void FORTRAN_NAME( FINDALLNEGATIVES ,findallnegatives )(
      CHFp_INT(count)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(interiorBox) );

#define FORT_FINDALLNEGATIVES FORTRAN_NAME( inlineFINDALLNEGATIVES, inlineFINDALLNEGATIVES)
#define FORTNT_FINDALLNEGATIVES FORTRAN_NAME( FINDALLNEGATIVES, findallnegatives)

inline void FORTRAN_NAME(inlineFINDALLNEGATIVES, inlineFINDALLNEGATIVES)(
      CHFp_INT(count)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(interiorBox) )
{
 CH_TIMELEAF("FORT_FINDALLNEGATIVES");
 FORTRAN_NAME( FINDALLNEGATIVES ,findallnegatives )(
      CHFt_INT(count)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(interiorBox) );
}
#endif  // GUARDFINDALLNEGATIVES 

}

#endif
