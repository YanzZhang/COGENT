/* Copyright (C) 1991-2015 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 7.0.0.  Version 7.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2012, plus Amendments 1 (published
   on April, 2013) and 2 (not yet published as of February, 2015).
   Additionally, it includes the accelerated publication of U+20BD
   RUBLE SIGN.  Therefore Unicode 7.0.0 is between 10646:2012 and
   10646:2014, and so we use the date ISO/IEC 10646:2012 Amd.1 was
   published.  */
/* We do not support C11 <threads.h>.  */
      subroutine COMPUTE_VEL_CELL_4D(
     & velCell
     & ,ivelCelllo0,ivelCelllo1,ivelCelllo2,ivelCelllo3
     & ,ivelCellhi0,ivelCellhi1,ivelCellhi2,ivelCellhi3
     & ,nvelCellcomp
     & ,velFace
     & ,ivelFacelo0,ivelFacelo1,ivelFacelo2,ivelFacelo3
     & ,ivelFacehi0,ivelFacehi1,ivelFacehi2,ivelFacehi3
     & ,nvelFacecomp
     & ,dir
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nvelCellcomp
      integer ivelCelllo0,ivelCelllo1,ivelCelllo2,ivelCelllo3
      integer ivelCellhi0,ivelCellhi1,ivelCellhi2,ivelCellhi3
      REAL*8 velCell(
     & ivelCelllo0:ivelCellhi0,
     & ivelCelllo1:ivelCellhi1,
     & ivelCelllo2:ivelCellhi2,
     & ivelCelllo3:ivelCellhi3,
     & 0:nvelCellcomp-1)
      integer nvelFacecomp
      integer ivelFacelo0,ivelFacelo1,ivelFacelo2,ivelFacelo3
      integer ivelFacehi0,ivelFacehi1,ivelFacehi2,ivelFacehi3
      REAL*8 velFace(
     & ivelFacelo0:ivelFacehi0,
     & ivelFacelo1:ivelFacehi1,
     & ivelFacelo2:ivelFacehi2,
     & ivelFacelo3:ivelFacehi3,
     & 0:nvelFacecomp-1)
      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, ii,jj,kk,ll, cfg_dim, comp
      double precision fac
      cfg_dim = 2
      fac = 1.0/4.0
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      kk = CHF_ID(2,dir)
      ll = CHF_ID(3,dir)
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         do comp = 0, cfg_dim - 1
           velCell(i,j,k,l,comp) = velCell(i,j,k,l,comp)
     & + fac * velFace(i,j,k,l,comp)
     & + fac * velFace(i+ ii,j+ jj,k+ kk,l+ ll,comp)
         enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_PERP_VEL_4D(
     & iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,result
     & ,iresultlo0,iresultlo1,iresultlo2,iresultlo3
     & ,iresulthi0,iresulthi1,iresulthi2,iresulthi3
     & ,nresultcomp
     & ,dfn
     & ,idfnlo0,idfnlo1,idfnlo2,idfnlo3
     & ,idfnhi0,idfnhi1,idfnhi2,idfnhi3
     & ,gkVel
     & ,igkVello0,igkVello1,igkVello2,igkVello3
     & ,igkVelhi0,igkVelhi1,igkVelhi2,igkVelhi3
     & ,ngkVelcomp
     & ,velCoords
     & ,ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
     & ,ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
     & ,nvelCoordscomp
     & ,b
     & ,iblo0,iblo1,iblo2,iblo3
     & ,ibhi0,ibhi1,ibhi2,ibhi3
     & ,nbcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nresultcomp
      integer iresultlo0,iresultlo1,iresultlo2,iresultlo3
      integer iresulthi0,iresulthi1,iresulthi2,iresulthi3
      REAL*8 result(
     & iresultlo0:iresulthi0,
     & iresultlo1:iresulthi1,
     & iresultlo2:iresulthi2,
     & iresultlo3:iresulthi3,
     & 0:nresultcomp-1)
      integer idfnlo0,idfnlo1,idfnlo2,idfnlo3
      integer idfnhi0,idfnhi1,idfnhi2,idfnhi3
      REAL*8 dfn(
     & idfnlo0:idfnhi0,
     & idfnlo1:idfnhi1,
     & idfnlo2:idfnhi2,
     & idfnlo3:idfnhi3)
      integer ngkVelcomp
      integer igkVello0,igkVello1,igkVello2,igkVello3
      integer igkVelhi0,igkVelhi1,igkVelhi2,igkVelhi3
      REAL*8 gkVel(
     & igkVello0:igkVelhi0,
     & igkVello1:igkVelhi1,
     & igkVello2:igkVelhi2,
     & igkVello3:igkVelhi3,
     & 0:ngkVelcomp-1)
      integer nvelCoordscomp
      integer ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
      integer ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
      REAL*8 velCoords(
     & ivelCoordslo0:ivelCoordshi0,
     & ivelCoordslo1:ivelCoordshi1,
     & ivelCoordslo2:ivelCoordshi2,
     & ivelCoordslo3:ivelCoordshi3,
     & 0:nvelCoordscomp-1)
      integer nbcomp
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL*8 b(
     & iblo0:ibhi0,
     & iblo1:ibhi1,
     & iblo2:ibhi2,
     & iblo3:ibhi3,
     & 0:nbcomp-1)
      integer i,j,k,l
      double precision bR, bphi, bZ, vpar_R, vpar_Phi, vpar_Z, vperp_R, 
     &vperp_Phi, vperp_Z
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          bR = b(i,j,iblo2,iblo3,0)
          bZ = b(i,j,iblo2,iblo3,2)
          vpar_R = velCoords(i,j,k,l,0) * bR
          vpar_Z = velCoords(i,j,k,l,0) * bZ
          vperp_R = gkVel(i,j,k,l,0) - vpar_R
          vperp_Z = gkVel(i,j,k,l,1) - vpar_Z
          result(i,j,k,l,0) = vperp_R * dfn(i,j,k,l)
          result(i,j,k,l,1) = vperp_Z * dfn(i,j,k,l)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_PAR_MOM_4D(
     & result
     & ,iresultlo0,iresultlo1,iresultlo2,iresultlo3
     & ,iresulthi0,iresulthi1,iresulthi2,iresulthi3
     & ,nresultcomp
     & ,parVel
     & ,iparVello0,iparVello1,iparVello2,iparVello3
     & ,iparVelhi0,iparVelhi1,iparVelhi2,iparVelhi3
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nresultcomp
      integer iresultlo0,iresultlo1,iresultlo2,iresultlo3
      integer iresulthi0,iresulthi1,iresulthi2,iresulthi3
      REAL*8 result(
     & iresultlo0:iresulthi0,
     & iresultlo1:iresulthi1,
     & iresultlo2:iresulthi2,
     & iresultlo3:iresulthi3,
     & 0:nresultcomp-1)
      integer iparVello0,iparVello1,iparVello2,iparVello3
      integer iparVelhi0,iparVelhi1,iparVelhi2,iparVelhi3
      REAL*8 parVel(
     & iparVello0:iparVelhi0,
     & iparVello1:iparVelhi1,
     & iparVello2:iparVelhi2,
     & iparVello3:iparVelhi3)
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, comp
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         do comp = 0, nresultcomp-1
           result(i,j,k,l,comp) = result(i,j,k,l,comp) * parVel(i,j,k,l)
         enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_PRESSURE_4D(
     & result
     & ,iresultlo0,iresultlo1,iresultlo2,iresultlo3
     & ,iresulthi0,iresulthi1,iresulthi2,iresulthi3
     & ,nresultcomp
     & ,velCoords
     & ,ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
     & ,ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
     & ,nvelCoordscomp
     & ,vparShift
     & ,ivparShiftlo0,ivparShiftlo1,ivparShiftlo2,ivparShiftlo3
     & ,ivparShifthi0,ivparShifthi1,ivparShifthi2,ivparShifthi3
     & ,B
     & ,iBlo0,iBlo1,iBlo2,iBlo3
     & ,iBhi0,iBhi1,iBhi2,iBhi3
     & ,mass
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nresultcomp
      integer iresultlo0,iresultlo1,iresultlo2,iresultlo3
      integer iresulthi0,iresulthi1,iresulthi2,iresulthi3
      REAL*8 result(
     & iresultlo0:iresulthi0,
     & iresultlo1:iresulthi1,
     & iresultlo2:iresulthi2,
     & iresultlo3:iresulthi3,
     & 0:nresultcomp-1)
      integer nvelCoordscomp
      integer ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
      integer ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
      REAL*8 velCoords(
     & ivelCoordslo0:ivelCoordshi0,
     & ivelCoordslo1:ivelCoordshi1,
     & ivelCoordslo2:ivelCoordshi2,
     & ivelCoordslo3:ivelCoordshi3,
     & 0:nvelCoordscomp-1)
      integer ivparShiftlo0,ivparShiftlo1,ivparShiftlo2,ivparShiftlo3
      integer ivparShifthi0,ivparShifthi1,ivparShifthi2,ivparShifthi3
      REAL*8 vparShift(
     & ivparShiftlo0:ivparShifthi0,
     & ivparShiftlo1:ivparShifthi1,
     & ivparShiftlo2:ivparShifthi2,
     & ivparShiftlo3:ivparShifthi3)
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL*8 B(
     & iBlo0:iBhi0,
     & iBlo1:iBhi1,
     & iBlo2:iBhi2,
     & iBlo3:iBhi3)
      REAL*8 mass
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, comp
      double precision vperp2, vparLoc, v2
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         vperp2 = velCoords(i,j,k,l,1) * B(i,j,iblo2,iblo3)
         vparLoc = velCoords(i,j,k,l,0) - vparShift(i,j,ivparShiftlo2,iv
     &parShiftlo3)
         v2 = mass * vparLoc**2 + vperp2
         do comp = 0, nresultcomp-1
           result(i,j,k,l,comp) = result(i,j,k,l,comp) * v2 / 3.0
         enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_MAXWELLIAN_4D(
     & result
     & ,iresultlo0,iresultlo1,iresultlo2,iresultlo3
     & ,iresulthi0,iresulthi1,iresulthi2,iresulthi3
     & ,nresultcomp
     & ,velCoords
     & ,ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
     & ,ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
     & ,nvelCoordscomp
     & ,vparShift
     & ,ivparShiftlo0,ivparShiftlo1,ivparShiftlo2,ivparShiftlo3
     & ,ivparShifthi0,ivparShifthi1,ivparShifthi2,ivparShifthi3
     & ,n
     & ,inlo0,inlo1,inlo2,inlo3
     & ,inhi0,inhi1,inhi2,inhi3
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & ,B
     & ,iBlo0,iBlo1,iBlo2,iBlo3
     & ,iBhi0,iBhi1,iBhi2,iBhi3
     & ,mass
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nresultcomp
      integer iresultlo0,iresultlo1,iresultlo2,iresultlo3
      integer iresulthi0,iresulthi1,iresulthi2,iresulthi3
      REAL*8 result(
     & iresultlo0:iresulthi0,
     & iresultlo1:iresulthi1,
     & iresultlo2:iresulthi2,
     & iresultlo3:iresulthi3,
     & 0:nresultcomp-1)
      integer nvelCoordscomp
      integer ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
      integer ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
      REAL*8 velCoords(
     & ivelCoordslo0:ivelCoordshi0,
     & ivelCoordslo1:ivelCoordshi1,
     & ivelCoordslo2:ivelCoordshi2,
     & ivelCoordslo3:ivelCoordshi3,
     & 0:nvelCoordscomp-1)
      integer ivparShiftlo0,ivparShiftlo1,ivparShiftlo2,ivparShiftlo3
      integer ivparShifthi0,ivparShifthi1,ivparShifthi2,ivparShifthi3
      REAL*8 vparShift(
     & ivparShiftlo0:ivparShifthi0,
     & ivparShiftlo1:ivparShifthi1,
     & ivparShiftlo2:ivparShifthi2,
     & ivparShiftlo3:ivparShifthi3)
      integer inlo0,inlo1,inlo2,inlo3
      integer inhi0,inhi1,inhi2,inhi3
      REAL*8 n(
     & inlo0:inhi0,
     & inlo1:inhi1,
     & inlo2:inhi2,
     & inlo3:inhi3)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL*8 T(
     & iTlo0:iThi0,
     & iTlo1:iThi1,
     & iTlo2:iThi2,
     & iTlo3:iThi3)
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL*8 B(
     & iBlo0:iBhi0,
     & iBlo1:iBhi1,
     & iBlo2:iBhi2,
     & iBlo3:iBhi3)
      REAL*8 mass
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, comp
      double precision v_parallel, mu, eparnorm, munorm, val, factor
      factor = sqrt((3.14159265358979323846264338327950288D0) * ((2.0d0)
     &/mass)**(3.0d0));
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         v_parallel = velCoords(i,j,k,l,0)
         mu = velCoords(i,j,k,l,1)
         eparnorm = ((1.0d0)/(2.0d0)) * mass * (v_parallel-vparShift(i,j
     &,ivparShiftlo2,ivparShiftlo3))**(2.0d0) / T(i,j,iTlo2,iTlo3)
         munorm = ((1.0d0)/(2.0d0)) * B(i,j,iblo2,iblo3) * mu / T(i,j,iT
     &lo2,iTlo3)
         val = exp( -( eparnorm + munorm ) )
         val = val * n(i,j,inlo2,inlo3) / ( factor * T(i,j,iTlo2,iTlo3)*
     &*((3.0d0)/(2.0d0)) )
         do comp = 0, nresultcomp-1
           result(i,j,k,l,comp) = val
         enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
