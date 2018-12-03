#ifndef _MAPPEDLIMITERF_F_H_
#define _MAPPEDLIMITERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCENTEREDLAPLACIAN 
#define GUARDCENTEREDLAPLACIAN 
// Prototype for Fortran procedure CENTEREDLAPLACIAN ...
//
void FORTRAN_NAME( CENTEREDLAPLACIAN ,centeredlaplacian )(
      CHFp_FRA(lapPhi)
      ,CHFp_CONST_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(dx) );

#define FORT_CENTEREDLAPLACIAN FORTRAN_NAME( inlineCENTEREDLAPLACIAN, inlineCENTEREDLAPLACIAN)
#define FORTNT_CENTEREDLAPLACIAN FORTRAN_NAME( CENTEREDLAPLACIAN, centeredlaplacian)

inline void FORTRAN_NAME(inlineCENTEREDLAPLACIAN, inlineCENTEREDLAPLACIAN)(
      CHFp_FRA(lapPhi)
      ,CHFp_CONST_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_CENTEREDLAPLACIAN");
 FORTRAN_NAME( CENTEREDLAPLACIAN ,centeredlaplacian )(
      CHFt_FRA(lapPhi)
      ,CHFt_CONST_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDCENTEREDLAPLACIAN 

#ifndef GUARDCCLAPLACIAN 
#define GUARDCCLAPLACIAN 
// Prototype for Fortran procedure CCLAPLACIAN ...
//
void FORTRAN_NAME( CCLAPLACIAN ,cclaplacian )(
      CHFp_FRA(lapPhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(dx) );

#define FORT_CCLAPLACIAN FORTRAN_NAME( inlineCCLAPLACIAN, inlineCCLAPLACIAN)
#define FORTNT_CCLAPLACIAN FORTRAN_NAME( CCLAPLACIAN, cclaplacian)

inline void FORTRAN_NAME(inlineCCLAPLACIAN, inlineCCLAPLACIAN)(
      CHFp_FRA(lapPhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_CCLAPLACIAN");
 FORTRAN_NAME( CCLAPLACIAN ,cclaplacian )(
      CHFt_FRA(lapPhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDCCLAPLACIAN 

#ifndef GUARDLIMITFACEVALUES 
#define GUARDLIMITFACEVALUES 
// Prototype for Fortran procedure LIMITFACEVALUES ...
//
void FORTRAN_NAME( LIMITFACEVALUES ,limitfacevalues )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA(centeredLap)
      ,CHFp_CONST_FRA(LRLap)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(limitC) );

#define FORT_LIMITFACEVALUES FORTRAN_NAME( inlineLIMITFACEVALUES, inlineLIMITFACEVALUES)
#define FORTNT_LIMITFACEVALUES FORTRAN_NAME( LIMITFACEVALUES, limitfacevalues)

inline void FORTRAN_NAME(inlineLIMITFACEVALUES, inlineLIMITFACEVALUES)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA(centeredLap)
      ,CHFp_CONST_FRA(LRLap)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(limitC) )
{
 CH_TIMELEAF("FORT_LIMITFACEVALUES");
 FORTRAN_NAME( LIMITFACEVALUES ,limitfacevalues )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA(centeredLap)
      ,CHFt_CONST_FRA(LRLap)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(limitC) );
}
#endif  // GUARDLIMITFACEVALUES 

#ifndef GUARDCOMPUTEA6 
#define GUARDCOMPUTEA6 
// Prototype for Fortran procedure COMPUTEA6 ...
//
void FORTRAN_NAME( COMPUTEA6 ,computea6 )(
      CHFp_FRA(a6)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA(facePhi)
      ,CHFp_BOX(ccbox)
      ,CHFp_CONST_INT(dir) );

#define FORT_COMPUTEA6 FORTRAN_NAME( inlineCOMPUTEA6, inlineCOMPUTEA6)
#define FORTNT_COMPUTEA6 FORTRAN_NAME( COMPUTEA6, computea6)

inline void FORTRAN_NAME(inlineCOMPUTEA6, inlineCOMPUTEA6)(
      CHFp_FRA(a6)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA(facePhi)
      ,CHFp_BOX(ccbox)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_COMPUTEA6");
 FORTRAN_NAME( COMPUTEA6 ,computea6 )(
      CHFt_FRA(a6)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA(facePhi)
      ,CHFt_BOX(ccbox)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDCOMPUTEA6 

#ifndef GUARDLIMITD2A 
#define GUARDLIMITD2A 
// Prototype for Fortran procedure LIMITD2a ...
//
void FORTRAN_NAME( LIMITD2A ,limitd2a )(
      CHFp_FRA(D2PhiLim)
      ,CHFp_CONST_FRA(D2Phi)
      ,CHFp_CONST_FRA(lapPhi)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(limitC) );

#define FORT_LIMITD2A FORTRAN_NAME( inlineLIMITD2A, inlineLIMITD2A)
#define FORTNT_LIMITD2A FORTRAN_NAME( LIMITD2A, limitd2a)

inline void FORTRAN_NAME(inlineLIMITD2A, inlineLIMITD2A)(
      CHFp_FRA(D2PhiLim)
      ,CHFp_CONST_FRA(D2Phi)
      ,CHFp_CONST_FRA(lapPhi)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(limitC) )
{
 CH_TIMELEAF("FORT_LIMITD2A");
 FORTRAN_NAME( LIMITD2A ,limitd2a )(
      CHFt_FRA(D2PhiLim)
      ,CHFt_CONST_FRA(D2Phi)
      ,CHFt_CONST_FRA(lapPhi)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(limitC) );
}
#endif  // GUARDLIMITD2A 

#ifndef GUARDLEFTRIGHTSTATES 
#define GUARDLEFTRIGHTSTATES 
// Prototype for Fortran procedure LEFTRIGHTSTATES ...
//
void FORTRAN_NAME( LEFTRIGHTSTATES ,leftrightstates )(
      CHFp_FRA(phiLeft)
      ,CHFp_FRA(phiRight)
      ,CHFp_CONST_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA(D2PhiLim)
      ,CHFp_CONST_FRA(D2Phi)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir) );

#define FORT_LEFTRIGHTSTATES FORTRAN_NAME( inlineLEFTRIGHTSTATES, inlineLEFTRIGHTSTATES)
#define FORTNT_LEFTRIGHTSTATES FORTRAN_NAME( LEFTRIGHTSTATES, leftrightstates)

inline void FORTRAN_NAME(inlineLEFTRIGHTSTATES, inlineLEFTRIGHTSTATES)(
      CHFp_FRA(phiLeft)
      ,CHFp_FRA(phiRight)
      ,CHFp_CONST_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA(D2PhiLim)
      ,CHFp_CONST_FRA(D2Phi)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_LEFTRIGHTSTATES");
 FORTRAN_NAME( LEFTRIGHTSTATES ,leftrightstates )(
      CHFt_FRA(phiLeft)
      ,CHFt_FRA(phiRight)
      ,CHFt_CONST_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA(D2PhiLim)
      ,CHFt_CONST_FRA(D2Phi)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDLEFTRIGHTSTATES 

#ifndef GUARDSELECTUPWIND 
#define GUARDSELECTUPWIND 
// Prototype for Fortran procedure SELECTUPWIND ...
//
void FORTRAN_NAME( SELECTUPWIND ,selectupwind )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(phiLeft)
      ,CHFp_CONST_FRA(phiRight)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox) );

#define FORT_SELECTUPWIND FORTRAN_NAME( inlineSELECTUPWIND, inlineSELECTUPWIND)
#define FORTNT_SELECTUPWIND FORTRAN_NAME( SELECTUPWIND, selectupwind)

inline void FORTRAN_NAME(inlineSELECTUPWIND, inlineSELECTUPWIND)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(phiLeft)
      ,CHFp_CONST_FRA(phiRight)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox) )
{
 CH_TIMELEAF("FORT_SELECTUPWIND");
 FORTRAN_NAME( SELECTUPWIND ,selectupwind )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(phiLeft)
      ,CHFt_CONST_FRA(phiRight)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox) );
}
#endif  // GUARDSELECTUPWIND 

}

#endif
