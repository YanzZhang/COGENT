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
      SUBROUTINE EVAL_DV_4D(
     & Dv
     & ,iDvlo0,iDvlo1,iDvlo2,iDvlo3
     & ,iDvhi0,iDvhi1,iDvhi2,iDvhi3
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,nphicomp
     & ,bmag
     & ,ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
     & ,ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & ,mTP
     & ,mFP
     & ,Nvpar
     & ,Nmu
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iDvlo0,iDvlo1,iDvlo2,iDvlo3
      integer iDvhi0,iDvhi1,iDvhi2,iDvhi3
      REAL*8 Dv(
     & iDvlo0:iDvhi0,
     & iDvlo1:iDvhi1,
     & iDvlo2:iDvhi2,
     & iDvlo3:iDvhi3)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3,
     & 0:nphicomp-1)
      integer ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
      integer ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
      REAL*8 bmag(
     & ibmaglo0:ibmaghi0,
     & ibmaglo1:ibmaghi1,
     & ibmaglo2:ibmaghi2,
     & ibmaglo3:ibmaghi3)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      REAL*8 mTP
      REAL*8 mFP
      integer Nvpar
      integer Nmu
      INTEGER i,j,k,l
      DOUBLE PRECISION dv_inv, b, h_v
      DOUBLE PRECISION, PARAMETER :: c1 = 1.0_8 / 24.0_8,
     & c2 = 9.0_8 / 8.0_8
      dv_inv = (1.0_8/dx(0))
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        IF ((k .EQ. -Nvpar/2) .OR. (k .EQ. Nvpar/2)) THEN
          Dv(i,j,k,l) = 0.0_8
        ELSE
          h_v = dv_inv * (
     & c1 * phi(i,j,k-2,l,0)
     & - c2 * phi(i,j,k-1,l,0)
     & + c2 * phi(i,j,k ,l,0)
     & - c1 * phi(i,j,k+1,l,0)
     & )
          Dv(i,j,k,l) = (mTP/mFP) * h_v;
        ENDIF
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
      SUBROUTINE EVAL_DMU_4D(
     & Dmu
     & ,iDmulo0,iDmulo1,iDmulo2,iDmulo3
     & ,iDmuhi0,iDmuhi1,iDmuhi2,iDmuhi3
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,nphicomp
     & ,bmag
     & ,ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
     & ,ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & ,mTP
     & ,mFP
     & ,Nvpar
     & ,Nmu
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iDmulo0,iDmulo1,iDmulo2,iDmulo3
      integer iDmuhi0,iDmuhi1,iDmuhi2,iDmuhi3
      REAL*8 Dmu(
     & iDmulo0:iDmuhi0,
     & iDmulo1:iDmuhi1,
     & iDmulo2:iDmuhi2,
     & iDmulo3:iDmuhi3)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3,
     & 0:nphicomp-1)
      integer ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
      integer ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
      REAL*8 bmag(
     & ibmaglo0:ibmaghi0,
     & ibmaglo1:ibmaghi1,
     & ibmaglo2:ibmaghi2,
     & ibmaglo3:ibmaghi3)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      REAL*8 mTP
      REAL*8 mFP
      integer Nvpar
      integer Nmu
      INTEGER i,j,k,l
      DOUBLE PRECISION dmu_inv, b, mu, h_mu
      DOUBLE PRECISION, PARAMETER :: c1 = 1.0_8 / 24.0_8,
     & c2 = 9.0_8 / 8.0_8
      dmu_inv = (1.0_8/dx(1))
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        IF ((l .EQ. 0) .OR. (l .EQ. Nmu)) THEN
          Dmu(i,j,k,l) = 0.0_8
        ELSE
          mu = l*dx(1)
          b = bmag(i,j,ibmaglo2,ibmaglo3)
          h_mu = dmu_inv * (
     & c1 * phi(i,j,k,l-2,0)
     & - c2 * phi(i,j,k,l-1,0)
     & + c2 * phi(i,j,k,l ,0)
     & - c1 * phi(i,j,k,l+1,0)
     & )
          Dmu(i,j,k,l) = 2.0_8 * (mTP/mFP) * (mTP/b) * mu * h_mu
        ENDIF
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
      SUBROUTINE EVAL_DVV_4D(
     & Dvv
     & ,iDvvlo0,iDvvlo1,iDvvlo2,iDvvlo3
     & ,iDvvhi0,iDvvhi1,iDvvhi2,iDvvhi3
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,nphicomp
     & ,bmag
     & ,ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
     & ,ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & ,mTP
     & ,mFP
     & ,Nvpar
     & ,Nmu
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iDvvlo0,iDvvlo1,iDvvlo2,iDvvlo3
      integer iDvvhi0,iDvvhi1,iDvvhi2,iDvvhi3
      REAL*8 Dvv(
     & iDvvlo0:iDvvhi0,
     & iDvvlo1:iDvvhi1,
     & iDvvlo2:iDvvhi2,
     & iDvvlo3:iDvvhi3)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3,
     & 0:nphicomp-1)
      integer ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
      integer ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
      REAL*8 bmag(
     & ibmaglo0:ibmaghi0,
     & ibmaglo1:ibmaghi1,
     & ibmaglo2:ibmaghi2,
     & ibmaglo3:ibmaghi3)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      REAL*8 mTP
      REAL*8 mFP
      integer Nvpar
      integer Nmu
      INTEGER i,j,k,l
      DOUBLE PRECISION dvsq_inv, g_vv1, g_vv2
      DOUBLE PRECISION, PARAMETER :: d1 = 1.0_8/12.0_8,
     & d2 = 4.0_8/3.0_8,
     & d3 = 5.0_8/2.0_8
      dvsq_inv = (1.0_8/dx(0)) * (1.0_8/dx(0))
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        IF ((k .EQ. -Nvpar/2) .OR. (k .EQ. Nvpar/2)) THEN
          Dvv(i,j,k,l) = 0.0_8
        ELSE
          g_vv1 = dvsq_inv * (
     & - d1 * phi(i,j,k+1,l,1)
     & + d2 * phi(i,j,k ,l,1)
     & - d3 * phi(i,j,k-1,l,1)
     & + d2 * phi(i,j,k-2,l,1)
     & - d1 * phi(i,j,k-3,l,1)
     & )
          g_vv2 = dvsq_inv * (
     & - d1 * phi(i,j,k+2,l,1)
     & + d2 * phi(i,j,k+1,l,1)
     & - d3 * phi(i,j,k ,l,1)
     & + d2 * phi(i,j,k-1,l,1)
     & - d1 * phi(i,j,k-2,l,1)
     & )
          Dvv(i,j,k,l) = -0.5_8 * (g_vv1 + g_vv2)
        ENDIF
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
      SUBROUTINE EVAL_DMUMU_4D(
     & Dmumu
     & ,iDmumulo0,iDmumulo1,iDmumulo2,iDmumulo3
     & ,iDmumuhi0,iDmumuhi1,iDmumuhi2,iDmumuhi3
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,nphicomp
     & ,bmag
     & ,ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
     & ,ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & ,mTP
     & ,mFP
     & ,Nvpar
     & ,Nmu
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iDmumulo0,iDmumulo1,iDmumulo2,iDmumulo3
      integer iDmumuhi0,iDmumuhi1,iDmumuhi2,iDmumuhi3
      REAL*8 Dmumu(
     & iDmumulo0:iDmumuhi0,
     & iDmumulo1:iDmumuhi1,
     & iDmumulo2:iDmumuhi2,
     & iDmumulo3:iDmumuhi3)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3,
     & 0:nphicomp-1)
      integer ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
      integer ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
      REAL*8 bmag(
     & ibmaglo0:ibmaghi0,
     & ibmaglo1:ibmaghi1,
     & ibmaglo2:ibmaghi2,
     & ibmaglo3:ibmaghi3)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      REAL*8 mTP
      REAL*8 mFP
      integer Nvpar
      integer Nmu
      INTEGER i,j,k,l
      DOUBLE PRECISION dmu_inv, dmusq_inv, b, mu, g_mu, g_mumu1, g_mumu2
     &, g_mumu
      DOUBLE PRECISION, PARAMETER :: d1 = 1.0_8/12.0_8,
     & d2 = 4.0_8/3.0_8,
     & d3 = 5.0_8/2.0_8,
     & c1 = 1.0_8 / 24.0_8,
     & c2 = 9.0_8 / 8.0_8
      dmu_inv = (1.0_8/dx(1))
      dmusq_inv = dmu_inv * dmu_inv
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        IF ((l .EQ. 0) .OR. (l .EQ. Nmu)) THEN
          Dmumu(i,j,k,l) = 0.0_8
        ELSE
          mu = l*dx(1)
          b = bmag(i,j,ibmaglo2,ibmaglo3)
          g_mu = dmu_inv * (
     & c1 * phi(i,j,k,l-2,1)
     & - c2 * phi(i,j,k,l-1,1)
     & + c2 * phi(i,j,k,l ,1)
     & - c1 * phi(i,j,k,l+1,1)
     & )
          g_mumu1 = dmusq_inv * (
     & - d1 * phi(i,j,k,l-3,1)
     & + d2 * phi(i,j,k,l-2,1)
     & - d3 * phi(i,j,k,l-1,1)
     & + d2 * phi(i,j,k,l ,1)
     & - d1 * phi(i,j,k,l+1,1)
     & )
          g_mumu2 = dmusq_inv * (
     & - d1 * phi(i,j,k,l-2,1)
     & + d2 * phi(i,j,k,l-1,1)
     & - d3 * phi(i,j,k,l ,1)
     & + d2 * phi(i,j,k,l+1,1)
     & - d1 * phi(i,j,k,l+2,1)
     & )
          g_mumu = 0.5_8 * (g_mumu1 + g_mumu2)
          Dmumu(i,j,k,l) = -(2.0_8*mTP*mTP*mu/(b*b)) * (2.0_8*mu*g_mumu 
     &+ g_mu)
        ENDIF
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
      SUBROUTINE EVAL_DVMU_4D(
     & Dvmu
     & ,iDvmulo0,iDvmulo1,iDvmulo2,iDvmulo3
     & ,iDvmuhi0,iDvmuhi1,iDvmuhi2,iDvmuhi3
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,nphicomp
     & ,bmag
     & ,ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
     & ,ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & ,mTP
     & ,mFP
     & ,Nvpar
     & ,Nmu
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iDvmulo0,iDvmulo1,iDvmulo2,iDvmulo3
      integer iDvmuhi0,iDvmuhi1,iDvmuhi2,iDvmuhi3
      REAL*8 Dvmu(
     & iDvmulo0:iDvmuhi0,
     & iDvmulo1:iDvmuhi1,
     & iDvmulo2:iDvmuhi2,
     & iDvmulo3:iDvmuhi3)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3,
     & 0:nphicomp-1)
      integer ibmaglo0,ibmaglo1,ibmaglo2,ibmaglo3
      integer ibmaghi0,ibmaghi1,ibmaghi2,ibmaghi3
      REAL*8 bmag(
     & ibmaglo0:ibmaghi0,
     & ibmaglo1:ibmaghi1,
     & ibmaglo2:ibmaghi2,
     & ibmaglo3:ibmaghi3)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      REAL*8 mTP
      REAL*8 mFP
      integer Nvpar
      integer Nmu
      INTEGER i,j,k,l
      DOUBLE PRECISION dvpar_inv, dmu_inv, g_vmu, b, mu,
     & g_mu_km3, g_mu_km2, g_mu_km1, g_mu_k, g_mu_kp1, g_mu_kp2
      DOUBLE PRECISION, PARAMETER :: c1 = 1.0_8 / 24.0_8,
     & c2 = 9.0_8 / 8.0_8
      dvpar_inv = (1.0_8/dx(0))
      dmu_inv = (1.0_8/dx(1))
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        IF ((l .EQ. 0) .OR. (l .EQ. Nmu) .OR. (k .EQ. -Nvpar/2) .OR. (k 
     &.EQ. Nvpar/2)) THEN
          Dvmu(i,j,k,l) = 0.0_8
        ELSE
          mu = l*dx(1)
          b = bmag(i,j,ibmaglo2,ibmaglo3)
          g_mu_km2 = dmu_inv * (
     & c1 * phi(i,j,k-2,l-2,1)
     & - c2 * phi(i,j,k-2,l-1,1)
     & + c2 * phi(i,j,k-2,l ,1)
     & - c1 * phi(i,j,k-2,l+1,1)
     & )
          g_mu_km1 = dmu_inv * (
     & c1 * phi(i,j,k-1,l-2,1)
     & - c2 * phi(i,j,k-1,l-1,1)
     & + c2 * phi(i,j,k-1,l ,1)
     & - c1 * phi(i,j,k-1,l+1,1)
     & )
          g_mu_k = dmu_inv * (
     & c1 * phi(i,j,k,l-2,1)
     & - c2 * phi(i,j,k,l-1,1)
     & + c2 * phi(i,j,k,l ,1)
     & - c1 * phi(i,j,k,l+1,1)
     & )
          g_mu_kp1 = dmu_inv * (
     & c1 * phi(i,j,k+1,l-2,1)
     & - c2 * phi(i,j,k+1,l-1,1)
     & + c2 * phi(i,j,k+1,l ,1)
     & - c1 * phi(i,j,k+1,l+1,1)
     & )
          g_vmu = dvpar_inv * (
     & c1 * g_mu_km2
     & - c2 * g_mu_km1
     & + c2 * g_mu_k
     & - c1 * g_mu_kp1
     & )
          Dvmu(i,j,k,l) = -2.0_8 * (mTP*mu/b) * g_vmu
        ENDIF
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
      SUBROUTINE FLUX_VPAR_4D(
     & flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
     & ,ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
     & ,nfluxcomp
     & ,Dv
     & ,iDvlo0,iDvlo1,iDvlo2,iDvlo3
     & ,iDvhi0,iDvhi1,iDvhi2,iDvhi3
     & ,Dvv
     & ,iDvvlo0,iDvvlo1,iDvvlo2,iDvvlo3
     & ,iDvvhi0,iDvvhi1,iDvvhi2,iDvvhi3
     & ,Dvmu
     & ,iDvmulo0,iDvmulo1,iDvmulo2,iDvmulo3
     & ,iDvmuhi0,iDvmuhi1,iDvmuhi2,iDvmuhi3
     & ,dfn
     & ,idfnlo0,idfnlo1,idfnlo2,idfnlo3
     & ,idfnhi0,idfnhi1,idfnhi2,idfnhi3
     & ,ndfncomp
     & ,dfn_f
     & ,idfn_flo0,idfn_flo1,idfn_flo2,idfn_flo3
     & ,idfn_fhi0,idfn_fhi1,idfn_fhi2,idfn_fhi3
     & ,ndfn_fcomp
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & ,Nvpar
     & ,Nmu
     & ,flag
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
      integer ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & ifluxlo3:ifluxhi3,
     & 0:nfluxcomp-1)
      integer iDvlo0,iDvlo1,iDvlo2,iDvlo3
      integer iDvhi0,iDvhi1,iDvhi2,iDvhi3
      REAL*8 Dv(
     & iDvlo0:iDvhi0,
     & iDvlo1:iDvhi1,
     & iDvlo2:iDvhi2,
     & iDvlo3:iDvhi3)
      integer iDvvlo0,iDvvlo1,iDvvlo2,iDvvlo3
      integer iDvvhi0,iDvvhi1,iDvvhi2,iDvvhi3
      REAL*8 Dvv(
     & iDvvlo0:iDvvhi0,
     & iDvvlo1:iDvvhi1,
     & iDvvlo2:iDvvhi2,
     & iDvvlo3:iDvvhi3)
      integer iDvmulo0,iDvmulo1,iDvmulo2,iDvmulo3
      integer iDvmuhi0,iDvmuhi1,iDvmuhi2,iDvmuhi3
      REAL*8 Dvmu(
     & iDvmulo0:iDvmuhi0,
     & iDvmulo1:iDvmuhi1,
     & iDvmulo2:iDvmuhi2,
     & iDvmulo3:iDvmuhi3)
      integer ndfncomp
      integer idfnlo0,idfnlo1,idfnlo2,idfnlo3
      integer idfnhi0,idfnhi1,idfnhi2,idfnhi3
      REAL*8 dfn(
     & idfnlo0:idfnhi0,
     & idfnlo1:idfnhi1,
     & idfnlo2:idfnhi2,
     & idfnlo3:idfnhi3,
     & 0:ndfncomp-1)
      integer ndfn_fcomp
      integer idfn_flo0,idfn_flo1,idfn_flo2,idfn_flo3
      integer idfn_fhi0,idfn_fhi1,idfn_fhi2,idfn_fhi3
      REAL*8 dfn_f(
     & idfn_flo0:idfn_fhi0,
     & idfn_flo1:idfn_fhi1,
     & idfn_flo2:idfn_fhi2,
     & idfn_flo3:idfn_fhi3,
     & 0:ndfn_fcomp-1)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      integer Nvpar
      integer Nmu
      integer flag
      INTEGER i,j,k,l
      DOUBLE PRECISION dvpar_inv, dmu_inv, alpha, dfn_face, D_vv, f_v, D
     &_vmu_times_f_mu,
     & f_mu_1, f_mu_2, f_mu_11, f_mu_12, f_mu_21, f_mu_22
      DOUBLE PRECISION mc
      DOUBLE PRECISION, PARAMETER :: c1 = 1.0_8 / 24.0_8,
     & c2 = 9.0_8 / 8.0_8
      dvpar_inv = 1.0_8 / dx(0)
      dmu_inv = 1.0_8 / dx(1)
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        alpha = Dv(i,j,k,l)
        dfn_face = dfn_f(i,j,k,l,0)
        D_vv = Dvv(i,j,k,l)
        f_v = dvpar_inv * (
     & c1 * dfn(i,j,k-2,l,0)
     & - c2 * dfn(i,j,k-1,l,0)
     & + c2 * dfn(i,j,k ,l,0)
     & - c1 * dfn(i,j,k+1,l,0)
     & )
        IF (flag .EQ. 0) THEN
          f_mu_11 = dmu_inv * (
     & c1 * dfn(i,j,k-1,l-2,0)
     & - c2 * dfn(i,j,k-1,l-1,0)
     & + c2 * dfn(i,j,k-1,l ,0)
     & - c1 * dfn(i,j,k-1,l+1,0)
     & )
          f_mu_12 = dmu_inv * (
     & c1 * dfn(i,j,k,l-2,0)
     & - c2 * dfn(i,j,k,l-1,0)
     & + c2 * dfn(i,j,k,l ,0)
     & - c1 * dfn(i,j,k,l+1,0)
     & )
          f_mu_21 = dmu_inv * (
     & c1 * dfn(i,j,k-1,l-1,0)
     & - c2 * dfn(i,j,k-1,l ,0)
     & + c2 * dfn(i,j,k-1,l+1,0)
     & - c1 * dfn(i,j,k-1,l+2,0)
     & )
          f_mu_22 = dmu_inv * (
     & c1 * dfn(i,j,k,l-1,0)
     & - c2 * dfn(i,j,k,l ,0)
     & + c2 * dfn(i,j,k,l+1,0)
     & - c1 * dfn(i,j,k,l+2,0)
     & )
          f_mu_1 = 0.5_8 * (f_mu_11 + f_mu_12)
          f_mu_2 = 0.5_8 * (f_mu_21 + f_mu_22)
          D_vmu_times_f_mu = 0.5_8 * ( Dvmu(i,j,k,l) * f_mu_1 + Dvmu(i,j
     &,k,l+1) * f_mu_2 )
        ELSE
          D_vmu_times_f_mu = mc(
     & Dvmu(i,j,k,l)
     & * mc(
     & dmu_inv*(dfn(i,j,k-1,l,0)-dfn(i,j,k-1,l-1,0)),
     & dmu_inv*(dfn(i,j,k,l,0)-dfn(i,j,k,l-1,0))
     & ),
     & Dvmu(i,j,k,l+1)
     & * mc(
     & dmu_inv*(dfn(i,j,k-1,l+1,0)-dfn(i,j,k-1,l,0)),
     & dmu_inv*(dfn(i,j,k,l+1,0)-dfn(i,j,k,l,0))
     & )
     & )
        ENDIF
        flux(i,j,k,l,0) = alpha*dfn_face
        flux(i,j,k,l,1) = D_vv*f_v + 2.0_8*D_vmu_times_f_mu
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
      SUBROUTINE FLUX_MU_4D(
     & flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
     & ,ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
     & ,nfluxcomp
     & ,Dmu
     & ,iDmulo0,iDmulo1,iDmulo2,iDmulo3
     & ,iDmuhi0,iDmuhi1,iDmuhi2,iDmuhi3
     & ,Dmumu
     & ,iDmumulo0,iDmumulo1,iDmumulo2,iDmumulo3
     & ,iDmumuhi0,iDmumuhi1,iDmumuhi2,iDmumuhi3
     & ,Dvmu
     & ,iDvmulo0,iDvmulo1,iDvmulo2,iDvmulo3
     & ,iDvmuhi0,iDvmuhi1,iDvmuhi2,iDvmuhi3
     & ,dfn
     & ,idfnlo0,idfnlo1,idfnlo2,idfnlo3
     & ,idfnhi0,idfnhi1,idfnhi2,idfnhi3
     & ,ndfncomp
     & ,dfn_f
     & ,idfn_flo0,idfn_flo1,idfn_flo2,idfn_flo3
     & ,idfn_fhi0,idfn_fhi1,idfn_fhi2,idfn_fhi3
     & ,ndfn_fcomp
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & ,Nvpar
     & ,Nmu
     & ,flag
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
      integer ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & ifluxlo3:ifluxhi3,
     & 0:nfluxcomp-1)
      integer iDmulo0,iDmulo1,iDmulo2,iDmulo3
      integer iDmuhi0,iDmuhi1,iDmuhi2,iDmuhi3
      REAL*8 Dmu(
     & iDmulo0:iDmuhi0,
     & iDmulo1:iDmuhi1,
     & iDmulo2:iDmuhi2,
     & iDmulo3:iDmuhi3)
      integer iDmumulo0,iDmumulo1,iDmumulo2,iDmumulo3
      integer iDmumuhi0,iDmumuhi1,iDmumuhi2,iDmumuhi3
      REAL*8 Dmumu(
     & iDmumulo0:iDmumuhi0,
     & iDmumulo1:iDmumuhi1,
     & iDmumulo2:iDmumuhi2,
     & iDmumulo3:iDmumuhi3)
      integer iDvmulo0,iDvmulo1,iDvmulo2,iDvmulo3
      integer iDvmuhi0,iDvmuhi1,iDvmuhi2,iDvmuhi3
      REAL*8 Dvmu(
     & iDvmulo0:iDvmuhi0,
     & iDvmulo1:iDvmuhi1,
     & iDvmulo2:iDvmuhi2,
     & iDvmulo3:iDvmuhi3)
      integer ndfncomp
      integer idfnlo0,idfnlo1,idfnlo2,idfnlo3
      integer idfnhi0,idfnhi1,idfnhi2,idfnhi3
      REAL*8 dfn(
     & idfnlo0:idfnhi0,
     & idfnlo1:idfnhi1,
     & idfnlo2:idfnhi2,
     & idfnlo3:idfnhi3,
     & 0:ndfncomp-1)
      integer ndfn_fcomp
      integer idfn_flo0,idfn_flo1,idfn_flo2,idfn_flo3
      integer idfn_fhi0,idfn_fhi1,idfn_fhi2,idfn_fhi3
      REAL*8 dfn_f(
     & idfn_flo0:idfn_fhi0,
     & idfn_flo1:idfn_fhi1,
     & idfn_flo2:idfn_fhi2,
     & idfn_flo3:idfn_fhi3,
     & 0:ndfn_fcomp-1)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      integer Nvpar
      integer Nmu
      integer flag
      INTEGER i,j,k,l
      DOUBLE PRECISION dvpar_inv, dmu_inv, alpha, dfn_face, D_mumu, f_mu
     &, D_vmu_times_f_v,
     & f_v_1, f_v_2, f_v_11, f_v_12, f_v_21, f_v_22
      DOUBLE PRECISION mc
      DOUBLE PRECISION, PARAMETER :: c1 = 1.0_8 / 24.0_8,
     & c2 = 9.0_8 / 8.0_8
      dvpar_inv = 1.0_8 / dx(0)
      dmu_inv = 1.0_8 / dx(1)
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        alpha = Dmu(i,j,k,l)
        dfn_face = dfn_f(i,j,k,l,0)
        D_mumu = Dmumu(i,j,k,l)
        f_mu = dmu_inv * (
     & c1 * dfn(i,j,k,l-2,0)
     & - c2 * dfn(i,j,k,l-1,0)
     & + c2 * dfn(i,j,k,l ,0)
     & - c1 * dfn(i,j,k,l+1,0)
     & )
        IF (flag .EQ. 0) THEN
          f_v_11 = dvpar_inv * (
     & c1 * dfn(i,j,k-2,l-1,0)
     & - c2 * dfn(i,j,k-1,l-1,0)
     & + c2 * dfn(i,j,k ,l-1,0)
     & - c1 * dfn(i,j,k+1,l-1,0)
     & )
          f_v_12 = dvpar_inv * (
     & c1 * dfn(i,j,k-2,l,0)
     & - c2 * dfn(i,j,k-1,l,0)
     & + c2 * dfn(i,j,k ,l,0)
     & - c1 * dfn(i,j,k+1,l,0)
     & )
          f_v_21 = dvpar_inv * (
     & c1 * dfn(i,j,k-1,l-1,0)
     & - c2 * dfn(i,j,k ,l-1,0)
     & + c2 * dfn(i,j,k+1,l-1,0)
     & - c1 * dfn(i,j,k+2,l-1,0)
     & )
          f_v_22 = dvpar_inv * (
     & c1 * dfn(i,j,k-1,l,0)
     & - c2 * dfn(i,j,k ,l,0)
     & + c2 * dfn(i,j,k+1,l,0)
     & - c1 * dfn(i,j,k+2,l,0)
     & )
          f_v_1 = 0.5_8 * (f_v_11 + f_v_12)
          f_v_2 = 0.5_8 * (f_v_21 + f_v_22)
          D_vmu_times_f_v = 0.5_8 * ( Dvmu(i,j,k,l) * f_v_1 + Dvmu(i,j,k
     &+1,l) * f_v_2 )
        ELSE
          D_vmu_times_f_v = mc(
     & Dvmu(i,j,k,l)
     & * mc(
     & dvpar_inv * (dfn(i,j,k,l-1,0)-dfn(i,j,k-1,l-1,0)),
     & dvpar_inv * (dfn(i,j,k,l,0)-dfn(i,j,k-1,l,0))
     & ),
     & Dvmu(i,j,k+1,l)
     & * mc(
     & dvpar_inv * (dfn(i,j,k+1,l-1,0)-dfn(i,j,k,l-1,0)),
     & dvpar_inv * (dfn(i,j,k+1,l,0)-dfn(i,j,k,l,0))
     & )
     & )
        ENDIF
        flux(i,j,k,l,0) = 2.0_8 * alpha*dfn_face
        flux(i,j,k,l,1) = 2.0_8 * (D_vmu_times_f_v + 2.0_8*D_mumu*f_mu)
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
      SUBROUTINE EVAL_RHS_COMPONENTS_4D(
     & rhs
     & ,irhslo0,irhslo1,irhslo2,irhslo3
     & ,irhshi0,irhshi1,irhshi2,irhshi3
     & ,nrhscomp
     & ,flux_vpar
     & ,iflux_vparlo0,iflux_vparlo1,iflux_vparlo2,iflux_vparlo3
     & ,iflux_vparhi0,iflux_vparhi1,iflux_vparhi2,iflux_vparhi3
     & ,nflux_vparcomp
     & ,flux_mu
     & ,iflux_mulo0,iflux_mulo1,iflux_mulo2,iflux_mulo3
     & ,iflux_muhi0,iflux_muhi1,iflux_muhi2,iflux_muhi3
     & ,nflux_mucomp
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2,irhslo3
      integer irhshi0,irhshi1,irhshi2,irhshi3
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2,
     & irhslo3:irhshi3,
     & 0:nrhscomp-1)
      integer nflux_vparcomp
      integer iflux_vparlo0,iflux_vparlo1,iflux_vparlo2,iflux_vparlo3
      integer iflux_vparhi0,iflux_vparhi1,iflux_vparhi2,iflux_vparhi3
      REAL*8 flux_vpar(
     & iflux_vparlo0:iflux_vparhi0,
     & iflux_vparlo1:iflux_vparhi1,
     & iflux_vparlo2:iflux_vparhi2,
     & iflux_vparlo3:iflux_vparhi3,
     & 0:nflux_vparcomp-1)
      integer nflux_mucomp
      integer iflux_mulo0,iflux_mulo1,iflux_mulo2,iflux_mulo3
      integer iflux_muhi0,iflux_muhi1,iflux_muhi2,iflux_muhi3
      REAL*8 flux_mu(
     & iflux_mulo0:iflux_muhi0,
     & iflux_mulo1:iflux_muhi1,
     & iflux_mulo2:iflux_muhi2,
     & iflux_mulo3:iflux_muhi3,
     & 0:nflux_mucomp-1)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      INTEGER i,j,k,l
      DOUBLE PRECISION dvpar_inv, dmu_inv
      dvpar_inv = 1.0_8 / dx(0)
      dmu_inv = 1.0_8 / dx(1)
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        rhs(i,j,k,l,0) = dvpar_inv*(flux_vpar(i,j,k+1,l,0)-flux_vpar(i,j
     &,k,l,0))
     & + dmu_inv*(flux_mu(i,j,k,l+1,0)-flux_mu(i,j,k,l,0))
        rhs(i,j,k,l,1) = dvpar_inv*(flux_vpar(i,j,k+1,l,1)-flux_vpar(i,j
     &,k,l,1))
     & + dmu_inv*(flux_mu(i,j,k,l+1,1)-flux_mu(i,j,k,l,1))
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
      SUBROUTINE EVAL_RHS_4D(
     & rhs
     & ,irhslo0,irhslo1,irhslo2,irhslo3
     & ,irhshi0,irhshi1,irhshi2,irhshi3
     & ,flux_vpar
     & ,iflux_vparlo0,iflux_vparlo1,iflux_vparlo2,iflux_vparlo3
     & ,iflux_vparhi0,iflux_vparhi1,iflux_vparhi2,iflux_vparhi3
     & ,nflux_vparcomp
     & ,flux_mu
     & ,iflux_mulo0,iflux_mulo1,iflux_mulo2,iflux_mulo3
     & ,iflux_muhi0,iflux_muhi1,iflux_muhi2,iflux_muhi3
     & ,nflux_mucomp
     & ,e
     & ,ielo0,ielo1,ielo2,ielo3
     & ,iehi0,iehi1,iehi2,iehi3
     & ,igridlo0,igridlo1,igridlo2,igridlo3
     & ,igridhi0,igridhi1,igridhi2,igridhi3
     & ,dx
     & ,Nvpar
     & ,Nmu
     & ,nu
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irhslo0,irhslo1,irhslo2,irhslo3
      integer irhshi0,irhshi1,irhshi2,irhshi3
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2,
     & irhslo3:irhshi3)
      integer nflux_vparcomp
      integer iflux_vparlo0,iflux_vparlo1,iflux_vparlo2,iflux_vparlo3
      integer iflux_vparhi0,iflux_vparhi1,iflux_vparhi2,iflux_vparhi3
      REAL*8 flux_vpar(
     & iflux_vparlo0:iflux_vparhi0,
     & iflux_vparlo1:iflux_vparhi1,
     & iflux_vparlo2:iflux_vparhi2,
     & iflux_vparlo3:iflux_vparhi3,
     & 0:nflux_vparcomp-1)
      integer nflux_mucomp
      integer iflux_mulo0,iflux_mulo1,iflux_mulo2,iflux_mulo3
      integer iflux_muhi0,iflux_muhi1,iflux_muhi2,iflux_muhi3
      REAL*8 flux_mu(
     & iflux_mulo0:iflux_muhi0,
     & iflux_mulo1:iflux_muhi1,
     & iflux_mulo2:iflux_muhi2,
     & iflux_mulo3:iflux_muhi3,
     & 0:nflux_mucomp-1)
      integer ielo0,ielo1,ielo2,ielo3
      integer iehi0,iehi1,iehi2,iehi3
      REAL*8 e(
     & ielo0:iehi0,
     & ielo1:iehi1,
     & ielo2:iehi2,
     & ielo3:iehi3)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      REAL*8 dx(0:3)
      integer Nvpar
      integer Nmu
      REAL*8 nu
      INTEGER i,j,k,l
      DOUBLE PRECISION dvpar_inv, dmu_inv, efac
      dvpar_inv = 1.0_8 / dx(0)
      dmu_inv = 1.0_8 / dx(1)
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
        efac = e(i,j,ielo2,ielo3)
        rhs(i,j,k,l) = nu * dvpar_inv
     & * (flux_vpar(i,j,k+1,l,0) - flux_vpar(i,j,k,l,0))
     & + nu * dmu_inv
     & * (flux_mu(i,j,k,l+1,0) - flux_mu(i,j,k,l,0))
     & + nu * efac * dvpar_inv
     & * (flux_vpar(i,j,k+1,l,1) - flux_vpar(i,j,k,l,1))
     & + nu * efac * dmu_inv
     & * (flux_mu(i,j,k,l+1,1) - flux_mu(i,j,k,l,1))
      enddo
      enddo
      enddo
      enddo
      RETURN
      END SUBROUTINE
