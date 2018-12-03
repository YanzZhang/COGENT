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
      function NuD(x)
      implicit none
      double precision NuD, x, pi, G, F
      pi=3.14159265358979
      G=(erf(x)-x*2.0/sqrt(pi)*exp(-x*x))/(2.0*x*x)
      F=erf(x)
      NuD=1/(x*x*x)*(F-G)
      return
      end
      function NuS(x)
      implicit none
      double precision NuS, x, pi, G
      pi=3.14159265358979
      G=(erf(x)-x*2.0/sqrt(pi)*exp(-x*x))/(2.0*x*x)
      NuS=4.0*G/x
      return
      end
      function NuPar(x)
      implicit none
      double precision NuPar, x, pi, G
      pi=3.14159265358979
      G=(erf(x)-x*2.0/sqrt(pi)*exp(-x*x))/(2.0*x*x)
      NuPar=2.0*G/(x*x*x)
      return
      end
      subroutine EVALUATE_TP_LORENTZ_CONST_NUD_4D(
     & ClsFlux
     & ,iClsFluxlo0,iClsFluxlo1,iClsFluxlo2,iClsFluxlo3
     & ,iClsFluxhi0,iClsFluxhi1,iClsFluxhi2,iClsFluxhi3
     & ,nClsFluxcomp
     & ,f
     & ,iflo0,iflo1,iflo2,iflo3
     & ,ifhi0,ifhi1,ifhi2,ifhi3
     & ,b
     & ,iblo0,iblo1,iblo2,iblo3
     & ,ibhi0,ibhi1,ibhi2,ibhi3
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & ,dx
     & ,m
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nClsFluxcomp
      integer iClsFluxlo0,iClsFluxlo1,iClsFluxlo2,iClsFluxlo3
      integer iClsFluxhi0,iClsFluxhi1,iClsFluxhi2,iClsFluxhi3
      REAL*8 ClsFlux(
     & iClsFluxlo0:iClsFluxhi0,
     & iClsFluxlo1:iClsFluxhi1,
     & iClsFluxlo2:iClsFluxhi2,
     & iClsFluxlo3:iClsFluxhi3,
     & 0:nClsFluxcomp-1)
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & iflo2:ifhi2,
     & iflo3:ifhi3)
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL*8 b(
     & iblo0:ibhi0,
     & iblo1:ibhi1,
     & iblo2:ibhi2,
     & iblo3:ibhi3)
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL*8 dx(0:3)
      REAL*8 m
      integer i0,i1,i2,i3
      double precision fmu,fvpar
      double precision CfmuLorentz,CfvparLorentz
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
        fmu=(1.0/12.0/dx(1))*(8.0*(f(i0,i1,i2,i3+1)-f(i0,i1,i2,i3-1))-(f
     &(i0,i1,i2,i3+2)-f(i0,i1,i2,i3-2)))
        fvpar=(1.0/12.0/dx(0))*(8.0*(f(i0,i1,i2+1,i3)-f(i0,i1,i2-1,i3))-
     &(f(i0,i1,i2+2,i3)-f(i0,i1,i2-2,i3)))
        if ((i3.eq.0) .or. (i3.eq.1)) then
          fmu=1.0/dx(1)*(-25.0/12.0*f(i0,i1,i2,i3)+4.0*f(i0,i1,i2,i3+1)-
     &3.0*f(i0,i1,i2,i3+2)+4.0/3.0*f(i0,i1,i2,i3+3)-1.0/4.0*f(i0,i1,i2,i
     &3+4))
        endif
        CfvparLorentz = 0.5*(b(i0,i1,iblo2,iblo3)/m*dx(1)*(i3+0.5)*fvpar
     &-2.0*dx(0)*(i2+0.5)*dx(1)*(i3+0.5)*fmu)
        CfmuLorentz = 0.5*(4.0*m/b(i0,i1,iblo2,iblo3)*dx(0)*dx(0)*(i2+0.
     &5)*(i2+0.5)*(i3+0.5)*dx(1)*fmu-2.0*dx(0)*(i2+0.5)*dx(1)*(i3+0.5)*f
     &vpar)
        ClsFlux(i0,i1,i2,i3,0) = CfvparLorentz
        ClsFlux(i0,i1,i2,i3,1) = CfmuLorentz
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVALUATE_TP_LORENTZ_4D(
     & ClsFlux
     & ,iClsFluxlo0,iClsFluxlo1,iClsFluxlo2,iClsFluxlo3
     & ,iClsFluxhi0,iClsFluxhi1,iClsFluxhi2,iClsFluxhi3
     & ,nClsFluxcomp
     & ,f
     & ,iflo0,iflo1,iflo2,iflo3
     & ,ifhi0,ifhi1,ifhi2,ifhi3
     & ,bmag
     & ,ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
     & ,ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & ,dx
     & ,m
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nClsFluxcomp
      integer iClsFluxlo0,iClsFluxlo1,iClsFluxlo2,iClsFluxlo3
      integer iClsFluxhi0,iClsFluxhi1,iClsFluxhi2,iClsFluxhi3
      REAL*8 ClsFlux(
     & iClsFluxlo0:iClsFluxhi0,
     & iClsFluxlo1:iClsFluxhi1,
     & iClsFluxlo2:iClsFluxhi2,
     & iClsFluxlo3:iClsFluxhi3,
     & 0:nClsFluxcomp-1)
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & iflo2:ifhi2,
     & iflo3:ifhi3)
      integer ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
      integer ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
      REAL*8 bmag(
     & ibmaglo0:ibmaghi0,
     & ibmaglo1:ibmaghi1,
     & ibmaglo2:ibmaghi2,
     & ibmaglo3:ibmaghi3)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL*8 T(
     & iTlo0:iThi0,
     & iTlo1:iThi1,
     & iTlo2:iThi2,
     & iTlo3:iThi3)
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL*8 dx(0:3)
      REAL*8 m
      integer i0,i1,i2,i3
      double precision fmu,fvpar, nu_D, NuD, b , v_th, x
      double precision CfmuLorentz,CfvparLorentz
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
        b = bmag(i0,i1,ibmaglo2,ibmaglo3)
        v_th=sqrt(2.0*T(i0,i1,iTlo2,iTlo3)/m)
        x=sqrt((i2+0.5)*(i2+0.5)*dx(0)*dx(0)+(i3+0.5)*dx(1)*b/m)/v_th
        nu_D = NuD(x)
        fmu=(1.0/12.0/dx(1))*(8.0*(f(i0,i1,i2,i3+1)-f(i0,i1,i2,i3-1))-(f
     &(i0,i1,i2,i3+2)-f(i0,i1,i2,i3-2)))
        fvpar=(1.0/12.0/dx(0))*(8.0*(f(i0,i1,i2+1,i3)-f(i0,i1,i2-1,i3))-
     &(f(i0,i1,i2+2,i3)-f(i0,i1,i2-2,i3)))
        if ((i3.eq.0) .or. (i3.eq.1)) then
          fmu=1.0/dx(1)*(-25.0/12.0*f(i0,i1,i2,i3)+4.0*f(i0,i1,i2,i3+1)-
     &3.0*f(i0,i1,i2,i3+2)+4.0/3.0*f(i0,i1,i2,i3+3)-1.0/4.0*f(i0,i1,i2,i
     &3+4))
        endif
        CfvparLorentz = 0.5*(b/m*dx(1)*(i3+0.5)*nu_D*fvpar-2.0*dx(0)*(i2
     &+0.5)*dx(1)*(i3+0.5)*nu_D*fmu)
        CfmuLorentz = 0.5*(4.0*m/b*dx(0)*dx(0)*(i2+0.5)*(i2+0.5)*(i3+0.5
     &)*dx(1)*nu_D*fmu-2.0*dx(0)*(i2+0.5)*dx(1)*(i3+0.5)*nu_D*fvpar)
        ClsFlux(i0,i1,i2,i3,0) = CfvparLorentz
        ClsFlux(i0,i1,i2,i3,1) = CfmuLorentz
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVALUATE_TP_ENERG_DIFF_4D(
     & ClsFlux
     & ,iClsFluxlo0,iClsFluxlo1,iClsFluxlo2,iClsFluxlo3
     & ,iClsFluxhi0,iClsFluxhi1,iClsFluxhi2,iClsFluxhi3
     & ,nClsFluxcomp
     & ,f
     & ,iflo0,iflo1,iflo2,iflo3
     & ,ifhi0,ifhi1,ifhi2,ifhi3
     & ,bmag
     & ,ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
     & ,ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & ,dx
     & ,m
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nClsFluxcomp
      integer iClsFluxlo0,iClsFluxlo1,iClsFluxlo2,iClsFluxlo3
      integer iClsFluxhi0,iClsFluxhi1,iClsFluxhi2,iClsFluxhi3
      REAL*8 ClsFlux(
     & iClsFluxlo0:iClsFluxhi0,
     & iClsFluxlo1:iClsFluxhi1,
     & iClsFluxlo2:iClsFluxhi2,
     & iClsFluxlo3:iClsFluxhi3,
     & 0:nClsFluxcomp-1)
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & iflo2:ifhi2,
     & iflo3:ifhi3)
      integer ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
      integer ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
      REAL*8 bmag(
     & ibmaglo0:ibmaghi0,
     & ibmaglo1:ibmaghi1,
     & ibmaglo2:ibmaghi2,
     & ibmaglo3:ibmaghi3)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL*8 T(
     & iTlo0:iThi0,
     & iTlo1:iThi1,
     & iTlo2:iThi2,
     & iTlo3:iThi3)
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL*8 dx(0:3)
      REAL*8 m
      integer i0,i1,i2,i3
      double precision x, v_th, b, fmu,fvpar, nu_S, nu_Par, NuD, NuS, Nu
     &Par
      double precision CfmuEdiff,CfvparEdiff
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
        b = bmag(i0,i1,ibmaglo2,ibmaglo3)
        v_th=sqrt(2.0*T(i0,i1,iTlo2,iTlo3)/m)
        x=sqrt((i2+0.5)*(i2+0.5)*dx(0)*dx(0)+(i3+0.5)*dx(1)*b/m)/v_th
        nu_S=NuS(x)
        nu_Par=NuPar(x)
        fmu=(1.0/12.0/dx(1))*(8.0*(f(i0,i1,i2,i3+1)-f(i0,i1,i2,i3-1))-(f
     &(i0,i1,i2,i3+2)-f(i0,i1,i2,i3-2)))
        fvpar=(1.0/12.0/dx(0))*(8.0*(f(i0,i1,i2+1,i3)-f(i0,i1,i2-1,i3))-
     &(f(i0,i1,i2+2,i3)-f(i0,i1,i2-2,i3)))
        if ((i3.eq.0) .or. (i3.eq.1)) then
          fmu=1.0/dx(1)*(-25.0/12.0*f(i0,i1,i2,i3)+4.0*f(i0,i1,i2,i3+1)-
     &3.0*f(i0,i1,i2,i3+2)+4.0/3.0*f(i0,i1,i2,i3+3)-1.0/4.0*f(i0,i1,i2,i
     &3+4))
        endif
        ClsFlux(i0,i1,i2,i3,0) = ClsFlux(i0,i1,i2,i3,0) + 0.5*nu_S*f(i0,
     &i1,i2,i3)*(i2+0.5)*dx(0)+0.5*nu_Par*(i2+0.5)*dx(0)*(2.0*(i3+0.5)*d
     &x(1)*fmu+(i2+0.5)*dx(0)*fvpar)
        ClsFlux(i0,i1,i2,i3,1) = ClsFlux(i0,i1,i2,i3,1) + nu_S*f(i0,i1,i
     &2,i3)*(i3+0.5)*dx(1)+nu_Par*(i3+0.5)*dx(1)*(2.0*(i3+0.5)*dx(1)*fmu
     &+(i2+0.5)*dx(0)*fvpar)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVALUATE_COLL_FLUX_COMBINED_4D(
     & fluxRHS
     & ,ifluxRHSlo0,ifluxRHSlo1,ifluxRHSlo2,ifluxRHSlo3
     & ,ifluxRHShi0,ifluxRHShi1,ifluxRHShi2,ifluxRHShi3
     & ,nfluxRHScomp
     & ,dir
     & ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,fluxFace
     & ,ifluxFacelo0,ifluxFacelo1,ifluxFacelo2,ifluxFacelo3
     & ,ifluxFacehi0,ifluxFacehi1,ifluxFacehi2,ifluxFacehi3
     & ,nfluxFacecomp
     & ,fluxCell
     & ,ifluxCelllo0,ifluxCelllo1,ifluxCelllo2,ifluxCelllo3
     & ,ifluxCellhi0,ifluxCellhi1,ifluxCellhi2,ifluxCellhi3
     & ,nfluxCellcomp
     & ,Nvpar
     & ,Nmu
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfluxRHScomp
      integer ifluxRHSlo0,ifluxRHSlo1,ifluxRHSlo2,ifluxRHSlo3
      integer ifluxRHShi0,ifluxRHShi1,ifluxRHShi2,ifluxRHShi3
      REAL*8 fluxRHS(
     & ifluxRHSlo0:ifluxRHShi0,
     & ifluxRHSlo1:ifluxRHShi1,
     & ifluxRHSlo2:ifluxRHShi2,
     & ifluxRHSlo3:ifluxRHShi3,
     & 0:nfluxRHScomp-1)
      integer dir
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nfluxFacecomp
      integer ifluxFacelo0,ifluxFacelo1,ifluxFacelo2,ifluxFacelo3
      integer ifluxFacehi0,ifluxFacehi1,ifluxFacehi2,ifluxFacehi3
      REAL*8 fluxFace(
     & ifluxFacelo0:ifluxFacehi0,
     & ifluxFacelo1:ifluxFacehi1,
     & ifluxFacelo2:ifluxFacehi2,
     & ifluxFacelo3:ifluxFacehi3,
     & 0:nfluxFacecomp-1)
      integer nfluxCellcomp
      integer ifluxCelllo0,ifluxCelllo1,ifluxCelllo2,ifluxCelllo3
      integer ifluxCellhi0,ifluxCellhi1,ifluxCellhi2,ifluxCellhi3
      REAL*8 fluxCell(
     & ifluxCelllo0:ifluxCellhi0,
     & ifluxCelllo1:ifluxCellhi1,
     & ifluxCelllo2:ifluxCellhi2,
     & ifluxCelllo3:ifluxCellhi3,
     & 0:nfluxCellcomp-1)
      integer Nvpar
      integer Nmu
      integer i,j,k,l
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         fluxRHS(i,j,k,l,0) = (0.0d0)
         fluxRHS(i,j,k,l,1) = (0.0d0)
         fluxRHS(i,j,k,l,2) = fluxFace(i,j,k,l,0)
         fluxRHS(i,j,k,l,3) = fluxFace(i,j,k,l,1)
         if (l.eq.1) then
          fluxRHS(i,j,k,l,3) = 1.0/4.0 * fluxCell(i,j,k,l-1,1)
     & + 13.0/12.0 * fluxCell(i,j,k,l,1)
     & - 5.0/12.0 * fluxCell(i,j,k,l+1,1)
     & + 1.0/12.0 * fluxCell(i,j,k,l+2,1)
         endif
         if ((l.eq.0).or.(l.eq.Nmu)) then
          fluxRHS(i,j,k,l,3) = 0.0
         endif
         if ((k.eq.-Nvpar/2).or.(k.eq.Nvpar/2)) then
          fluxRHS(i,j,k,l,2) = 0.0
         endif
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVALUATE_FULL_ER_FLUX_4D(
     & dir
     & ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
     & ,ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
     & ,nfluxcomp
     & ,ERest
     & ,iERestlo0,iERestlo1,iERestlo2,iERestlo3
     & ,iEResthi0,iEResthi1,iEResthi2,iEResthi3
     & ,ENorm
     & ,iENormlo0,iENormlo1,iENormlo2,iENormlo3
     & ,iENormhi0,iENormhi1,iENormhi2,iENormhi3
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & ,b
     & ,iblo0,iblo1,iblo2,iblo3
     & ,ibhi0,ibhi1,ibhi2,ibhi3
     & ,dx
     & ,Nvpar
     & ,Nmu
     & ,m
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer dir
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
      integer ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & ifluxlo3:ifluxhi3,
     & 0:nfluxcomp-1)
      integer iERestlo0,iERestlo1,iERestlo2,iERestlo3
      integer iEResthi0,iEResthi1,iEResthi2,iEResthi3
      REAL*8 ERest(
     & iERestlo0:iEResthi0,
     & iERestlo1:iEResthi1,
     & iERestlo2:iEResthi2,
     & iERestlo3:iEResthi3)
      integer iENormlo0,iENormlo1,iENormlo2,iENormlo3
      integer iENormhi0,iENormhi1,iENormhi2,iENormhi3
      REAL*8 ENorm(
     & iENormlo0:iENormhi0,
     & iENormlo1:iENormhi1,
     & iENormlo2:iENormhi2,
     & iENormlo3:iENormhi3)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL*8 T(
     & iTlo0:iThi0,
     & iTlo1:iThi1,
     & iTlo2:iThi2,
     & iTlo3:iThi3)
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL*8 b(
     & iblo0:ibhi0,
     & iblo1:ibhi1,
     & iblo2:ibhi2,
     & iblo3:ibhi3)
      REAL*8 dx(0:3)
      integer Nvpar
      integer Nmu
      REAL*8 m
      integer i0,i1,i2,i3
      double precision x,v_th, Norm, nu_S, NuS
      do i3 = iboxlo3,iboxhi3
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
       flux(i0,i1,i2,i3,0) = (0.0d0)
       flux(i0,i1,i2,i3,1) = (0.0d0)
       Norm = ERest(i0,i1,iERestlo2,iERestlo3)/ENorm(i0,i1,iENormlo2,iEN
     &ormlo3)
       v_th=sqrt(2.0*T(i0,i1,iTlo2,iTlo3)/m)
       if ((i3.eq.0) .and. (i2.eq.0)) then
          flux(i0,i1,i2,i3,2) = 0.0
          flux(i0,i1,i2,i3,3) = 0.0
       else
         x=sqrt((i2)*(i2)*dx(0)*dx(0)+(i3+0.5)*dx(1)*b(i0,i1,iblo2,iblo3
     &)/m)/v_th
         nu_S=NuS(x)
         if ((i2.eq.-Nvpar/2).or.(i2.eq.Nvpar/2)) then
            flux(i0,i1,i2,i3,2) = 0.0
         else
           flux(i0,i1,i2,i3,2) = (i2)*dx(0)*0.5*v_th*v_th*nu_S*exp(-x*x)
     &*Norm
         endif
         x=sqrt((i2+0.5)*(i2+0.5)*dx(0)*dx(0)+(i3)*dx(1)*b(i0,i1,iblo2,i
     &blo3)/m)/v_th
         nu_S=NuS(x)
         if (i3.eq.Nmu) then
            flux(i0,i1,i2,i3,3) = 0.0
         else
           flux(i0,i1,i2,i3,3) = (i3)*dx(1)*v_th*v_th*nu_S*exp(-x*x)*Nor
     &m
         endif
       endif
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVALUATE_NORM_ER_FLUX_4D(
     & dir
     & ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
     & ,ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
     & ,nfluxcomp
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & ,b
     & ,iblo0,iblo1,iblo2,iblo3
     & ,ibhi0,ibhi1,ibhi2,ibhi3
     & ,dx
     & ,Nvpar
     & ,Nmu
     & ,m
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer dir
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
      integer ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & ifluxlo3:ifluxhi3,
     & 0:nfluxcomp-1)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL*8 T(
     & iTlo0:iThi0,
     & iTlo1:iThi1,
     & iTlo2:iThi2,
     & iTlo3:iThi3)
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL*8 b(
     & iblo0:ibhi0,
     & iblo1:ibhi1,
     & iblo2:ibhi2,
     & iblo3:ibhi3)
      REAL*8 dx(0:3)
      integer Nvpar
      integer Nmu
      REAL*8 m
      integer i0,i1,i2,i3
      double precision x,v_th, nu_S, NuS
      do i3 = iboxlo3,iboxhi3
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
       flux(i0,i1,i2,i3,0) = (0.0d0)
       flux(i0,i1,i2,i3,1) = (0.0d0)
       v_th=sqrt(2.0*T(i0,i1,iTlo2,iTlo3)/m)
       if ((i3.eq.0) .and. (i2.eq.0)) then
          flux(i0,i1,i2,i3,2) = 0.0
          flux(i0,i1,i2,i3,3) = 0.0
       else
         x=sqrt((i2)*(i2)*dx(0)*dx(0)+(i3+0.5)*dx(1)*b(i0,i1,iblo2,iblo3
     &)/m)/v_th
         nu_S=NuS(x)
         if ((i2.eq.-Nvpar/2).or.(i2.eq.Nvpar/2)) then
            flux(i0,i1,i2,i3,2) = 0.0
         else
            flux(i0,i1,i2,i3,2) = -i2*dx(0)*0.5*v_th*v_th*nu_S*exp(-x*x)
         endif
         x=sqrt((i2+0.5)*(i2+0.5)*dx(0)*dx(0)+(i3)*dx(1)*b(i0,i1,iblo2,i
     &blo3)/m)/v_th
         nu_S=NuS(x)
         if (i3.eq.Nmu) then
            flux(i0,i1,i2,i3,3) = 0.0
         else
            flux(i0,i1,i2,i3,3) = -i3*dx(1)*v_th*v_th*nu_S*exp(-x*x)
         endif
       endif
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_SC_CLS_FREQ_4D(
     & iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,cls_freq
     & ,icls_freqlo0,icls_freqlo1,icls_freqlo2,icls_freqlo3
     & ,icls_freqhi0,icls_freqhi1,icls_freqhi2,icls_freqhi3
     & ,n
     & ,inlo0,inlo1,inlo2,inlo3
     & ,inhi0,inhi1,inhi2,inhi3
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer icls_freqlo0,icls_freqlo1,icls_freqlo2,icls_freqlo3
      integer icls_freqhi0,icls_freqhi1,icls_freqhi2,icls_freqhi3
      REAL*8 cls_freq(
     & icls_freqlo0:icls_freqhi0,
     & icls_freqlo1:icls_freqhi1,
     & icls_freqlo2:icls_freqhi2,
     & icls_freqlo3:icls_freqhi3)
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
      integer i,j,k,l
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
       cls_freq(i,j,k,l) = n(i,j,inlo2,inlo3)
     & / (T(i,j,iTlo2,iTlo3)**(3.0/2.0))
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVAL_CONSDRAGDIFF_FLUX_4D(
     & iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,fluxes
     & ,ifluxeslo0,ifluxeslo1,ifluxeslo2,ifluxeslo3
     & ,ifluxeshi0,ifluxeshi1,ifluxeshi2,ifluxeshi3
     & ,nfluxescomp
     & ,fBJ
     & ,ifBJlo0,ifBJlo1,ifBJlo2,ifBJlo3
     & ,ifBJhi0,ifBJhi1,ifBJhi2,ifBJhi3
     & ,B
     & ,iBlo0,iBlo1,iBlo2,iBlo3
     & ,iBhi0,iBhi1,iBhi2,iBhi3
     & ,mass
     & ,dx
     & ,dir
     & ,Nvp
     & ,Nmu
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nfluxescomp
      integer ifluxeslo0,ifluxeslo1,ifluxeslo2,ifluxeslo3
      integer ifluxeshi0,ifluxeshi1,ifluxeshi2,ifluxeshi3
      REAL*8 fluxes(
     & ifluxeslo0:ifluxeshi0,
     & ifluxeslo1:ifluxeshi1,
     & ifluxeslo2:ifluxeshi2,
     & ifluxeslo3:ifluxeshi3,
     & 0:nfluxescomp-1)
      integer ifBJlo0,ifBJlo1,ifBJlo2,ifBJlo3
      integer ifBJhi0,ifBJhi1,ifBJhi2,ifBJhi3
      REAL*8 fBJ(
     & ifBJlo0:ifBJhi0,
     & ifBJlo1:ifBJhi1,
     & ifBJlo2:ifBJhi2,
     & ifBJlo3:ifBJhi3)
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL*8 B(
     & iBlo0:iBhi0,
     & iBlo1:iBhi1,
     & iBlo2:iBhi2,
     & iBlo3:iBhi3)
      REAL*8 mass
      REAL*8 dx(0:3)
      integer dir
      integer Nvp
      integer Nmu
      integer i,j,k,l, B2, B3, B4
      double precision fBJ_face, dfBJdvp, dfBJdmu
      B2 = iBlo2
      B3 = iBlo3
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
       if (dir==2) then
         if ((k==Nvp/2) .or. (k==-Nvp/2)) then
           fBJ_face = (0.0d0)
           dfBJdvp = (0.0d0)
         else
         if (k==Nvp/2-1) then
           fBJ_face = ( (3.0d0)*fBJ(i,j,k,l)+13.0d0*fBJ(i,j,k-1,l)-(5.0d
     &0)*fBJ(i,j,k-2,l)+fBJ(i,j,k-3,l) )/12.0d0
           dfBJdvp = ( (10.0d0)*fBJ(i,j,k,l)-(5.0d0)*fBJ(i,j,k-1,l)-(9.0
     &d0)*fBJ(i,j,k-2,l)+(5.0d0)*fBJ(i,j,k-3,l)-fBJ(i,j,k-4,l) )/12.0d0/
     &dx(2)
         else
         if (k==-Nvp/2+1) then
           fBJ_face = ( (3.0d0)*fBJ(i,j,k-1,l)+13.0d0*fBJ(i,j,k,l)-(5.0d
     &0)*fBJ(i,j,k+1,l)+fBJ(i,j,k+2,l) )/12.0d0
           dfBJdvp = (-(10.0d0)*fBJ(i,j,k-1,l)+(5.0d0)*fBJ(i,j,k,l)+(9.0
     &d0)*fBJ(i,j,k+1,l)-(5.0d0)*fBJ(i,j,k+2,l)+fBJ(i,j,k+3,l) )/12.0d0/
     &dx(2)
         else
           fBJ_face = ( (7.0d0)*(fBJ(i,j,k,l)+fBJ(i,j,k-1,l))-(fBJ(i,j,k
     &+1,l)+fBJ(i,j,k-2,l)) )/12.0d0
           dfBJdvp = (15.0d0*(fBJ(i,j,k,l)-fBJ(i,j,k-1,l))-(fBJ(i,j,k+1,
     &l)-fBJ(i,j,k-2,l)) )/12.0d0/dx(2)
         endif
         endif
         endif
         fluxes(i,j,k,l,0) = k*dx(2)*fBJ_face
         fluxes(i,j,k,l,1) = fBJ_face
         fluxes(i,j,k,l,2) = dfBJdvp/mass
       else
       if (dir==3) then
         if ((l==0) .or. (l==Nmu)) then
           fBJ_face = (0.0d0)
           dfBJdmu = (0.0d0)
         else
         if (l==Nmu-1) then
           fBJ_face = ( (3.0d0)*fBJ(i,j,k,l)+13.0d0*fBJ(i,j,k,l-1)-(5.0d
     &0)*fBJ(i,j,k,l-2)+fBJ(i,j,k,l-3) )/12.0d0
           dfBJdmu = ( (10.0d0)*fBJ(i,j,k,l)-(5.0d0)*fBJ(i,j,k,l-1)-(9.0
     &d0)*fBJ(i,j,k,l-2)+(5.0d0)*fBJ(i,j,k,l-3)-fBJ(i,j,k,l-4) )/12.0d0/
     &dx(3)
         else
         if (l==1) then
           fBJ_face = ( (3.0d0)*fBJ(i,j,k,l-1)+13.0d0*fBJ(i,j,k,l)-(5.0d
     &0)*fBJ(i,j,k,l+1)+fBJ(i,j,k,l+2) )/12.0d0
           dfBJdmu = (-(10.0d0)*fBJ(i,j,k,l-1)+(5.0d0)*fBJ(i,j,k,l)+(9.0
     &d0)*fBJ(i,j,k,l+1)-(5.0d0)*fBJ(i,j,k,l+2)+fBJ(i,j,k,l+3) )/12.0d0/
     &dx(3)
         else
           fBJ_face = ( (7.0d0)*(fBJ(i,j,k,l)+fBJ(i,j,k,l-1))-(fBJ(i,j,k
     &,l+1)+fBJ(i,j,k,l-2)) )/12.0d0
           dfBJdmu = (15.0d0*(fBJ(i,j,k,l)-fBJ(i,j,k,l-1))-(fBJ(i,j,k,l+
     &1)-fBJ(i,j,k,l-2)) )/12.0d0/dx(3)
         endif
         endif
         endif
         fluxes(i,j,k,l,0) = (2.0d0)*l*dx(3)*fBJ_face
         fluxes(i,j,k,l,1) = (0.0d0)
         fluxes(i,j,k,l,2) = (4.0d0)*l*dx(3)/B(i,j,B2,B3)*dfBJdmu
       else
         fluxes(i,j,k,l,0) = (0.0d0)
         fluxes(i,j,k,l,1) = (0.0d0)
         fluxes(i,j,k,l,2) = (0.0d0)
       endif
       endif
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVAL_CONS_UPAR_TEMP_4D(
     & iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,Upar
     & ,iUparlo0,iUparlo1,iUparlo2,iUparlo3
     & ,iUparhi0,iUparhi1,iUparhi2,iUparhi3
     & ,Temp
     & ,iTemplo0,iTemplo1,iTemplo2,iTemplo3
     & ,iTemphi0,iTemphi1,iTemphi2,iTemphi3
     & ,vmoms
     & ,ivmomslo0,ivmomslo1,ivmomslo2,ivmomslo3
     & ,ivmomshi0,ivmomshi1,ivmomshi2,ivmomshi3
     & ,nvmomscomp
     & ,pmoms
     & ,ipmomslo0,ipmomslo1,ipmomslo2,ipmomslo3
     & ,ipmomshi0,ipmomshi1,ipmomshi2,ipmomshi3
     & ,npmomscomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer iUparlo0,iUparlo1,iUparlo2,iUparlo3
      integer iUparhi0,iUparhi1,iUparhi2,iUparhi3
      REAL*8 Upar(
     & iUparlo0:iUparhi0,
     & iUparlo1:iUparhi1,
     & iUparlo2:iUparhi2,
     & iUparlo3:iUparhi3)
      integer iTemplo0,iTemplo1,iTemplo2,iTemplo3
      integer iTemphi0,iTemphi1,iTemphi2,iTemphi3
      REAL*8 Temp(
     & iTemplo0:iTemphi0,
     & iTemplo1:iTemphi1,
     & iTemplo2:iTemphi2,
     & iTemplo3:iTemphi3)
      integer nvmomscomp
      integer ivmomslo0,ivmomslo1,ivmomslo2,ivmomslo3
      integer ivmomshi0,ivmomshi1,ivmomshi2,ivmomshi3
      REAL*8 vmoms(
     & ivmomslo0:ivmomshi0,
     & ivmomslo1:ivmomshi1,
     & ivmomslo2:ivmomshi2,
     & ivmomslo3:ivmomshi3,
     & 0:nvmomscomp-1)
      integer npmomscomp
      integer ipmomslo0,ipmomslo1,ipmomslo2,ipmomslo3
      integer ipmomshi0,ipmomshi1,ipmomshi2,ipmomshi3
      REAL*8 pmoms(
     & ipmomslo0:ipmomshi0,
     & ipmomslo1:ipmomshi1,
     & ipmomslo2:ipmomshi2,
     & ipmomslo3:ipmomshi3,
     & 0:npmomscomp-1)
      integer i,j,k,l, L2, L3 ,L4
      double precision denom
      L2 = iUparlo2
      L3 = iUparlo3
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
       denom = pmoms(i,j,L2,L3,2)*vmoms(i,j,L2,L3,1)-pmoms(i,j,L2,L3,1)*
     &vmoms(i,j,L2,L3,2)
       Upar(i,j,L2,L3) = (pmoms(i,j,L2,L3,2)*vmoms(i,j,L2,L3,0)-pmoms(i,
     &j,L2,L3,0)*vmoms(i,j,L2,L3,2))/denom
       Temp(i,j,L2,L3) = (pmoms(i,j,L2,L3,1)*vmoms(i,j,L2,L3,0)-pmoms(i,
     &j,L2,L3,0)*vmoms(i,j,L2,L3,1))/denom
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVAL_CONSDRAGDIFF_RHS_4D(
     & iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,rhs_cls
     & ,irhs_clslo0,irhs_clslo1,irhs_clslo2,irhs_clslo3
     & ,irhs_clshi0,irhs_clshi1,irhs_clshi2,irhs_clshi3
     & ,nu
     & ,Upar
     & ,iUparlo0,iUparlo1,iUparlo2,iUparlo3
     & ,iUparhi0,iUparhi1,iUparhi2,iUparhi3
     & ,Temp
     & ,iTemplo0,iTemplo1,iTemplo2,iTemplo3
     & ,iTemphi0,iTemphi1,iTemphi2,iTemphi3
     & ,Jpsi
     & ,iJpsilo0,iJpsilo1,iJpsilo2,iJpsilo3
     & ,iJpsihi0,iJpsihi1,iJpsihi2,iJpsihi3
     & ,nJpsicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer irhs_clslo0,irhs_clslo1,irhs_clslo2,irhs_clslo3
      integer irhs_clshi0,irhs_clshi1,irhs_clshi2,irhs_clshi3
      REAL*8 rhs_cls(
     & irhs_clslo0:irhs_clshi0,
     & irhs_clslo1:irhs_clshi1,
     & irhs_clslo2:irhs_clshi2,
     & irhs_clslo3:irhs_clshi3)
      REAL*8 nu
      integer iUparlo0,iUparlo1,iUparlo2,iUparlo3
      integer iUparhi0,iUparhi1,iUparhi2,iUparhi3
      REAL*8 Upar(
     & iUparlo0:iUparhi0,
     & iUparlo1:iUparhi1,
     & iUparlo2:iUparhi2,
     & iUparlo3:iUparhi3)
      integer iTemplo0,iTemplo1,iTemplo2,iTemplo3
      integer iTemphi0,iTemphi1,iTemphi2,iTemphi3
      REAL*8 Temp(
     & iTemplo0:iTemphi0,
     & iTemplo1:iTemphi1,
     & iTemplo2:iTemphi2,
     & iTemplo3:iTemphi3)
      integer nJpsicomp
      integer iJpsilo0,iJpsilo1,iJpsilo2,iJpsilo3
      integer iJpsihi0,iJpsihi1,iJpsihi2,iJpsihi3
      REAL*8 Jpsi(
     & iJpsilo0:iJpsihi0,
     & iJpsilo1:iJpsihi1,
     & iJpsilo2:iJpsihi2,
     & iJpsilo3:iJpsihi3,
     & 0:nJpsicomp-1)
      integer i,j,k,l, L2, L3, L4
      L2 = iUparlo2
      L3 = iUparlo3
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
       rhs_cls(i,j,k,l) = nu*(Jpsi(i,j,k,l,0)
     & - Upar(i,j,L2,L3)*Jpsi(i,j,k,l,1)
     & + Temp(i,j,L2,L3)*Jpsi(i,j,k,l,2))
      enddo
      enddo
      enddo
      enddo
      return
      end
