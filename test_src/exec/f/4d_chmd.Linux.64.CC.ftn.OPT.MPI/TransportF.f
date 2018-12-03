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
      subroutine EVAL_ANOM_FLUX_4D(
     & dir
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & ,dx
     & ,lame_face
     & ,ilame_facelo0,ilame_facelo1,ilame_facelo2,ilame_facelo3
     & ,ilame_facehi0,ilame_facehi1,ilame_facehi2,ilame_facehi3
     & ,nlame_facecomp
     & ,D_cfg
     & ,iD_cfglo0,iD_cfglo1,iD_cfglo2,iD_cfglo3
     & ,iD_cfghi0,iD_cfghi1,iD_cfghi2,iD_cfghi3
     & ,nD_cfgcomp
     & ,mass
     & ,Nr
     & ,simpleDiffusion
     & ,fB
     & ,ifBlo0,ifBlo1,ifBlo2,ifBlo3
     & ,ifBhi0,ifBhi1,ifBhi2,ifBhi3
     & ,dfBdmu
     & ,idfBdmulo0,idfBdmulo1,idfBdmulo2,idfBdmulo3
     & ,idfBdmuhi0,idfBdmuhi1,idfBdmuhi2,idfBdmuhi3
     & ,B
     & ,iBlo0,iBlo1,iBlo2,iBlo3
     & ,iBhi0,iBhi1,iBhi2,iBhi3
     & ,N
     & ,iNlo0,iNlo1,iNlo2,iNlo3
     & ,iNhi0,iNhi1,iNhi2,iNhi3
     & ,U
     & ,iUlo0,iUlo1,iUlo2,iUlo3
     & ,iUhi0,iUhi1,iUhi2,iUhi3
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & ,C
     & ,iClo0,iClo1,iClo2,iClo3
     & ,iChi0,iChi1,iChi2,iChi3
     & ,P
     & ,iPlo0,iPlo1,iPlo2,iPlo3
     & ,iPhi0,iPhi1,iPhi2,iPhi3
     & ,flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
     & ,ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL*8 dx(0:3)
      integer nlame_facecomp
      integer ilame_facelo0,ilame_facelo1,ilame_facelo2,ilame_facelo3
      integer ilame_facehi0,ilame_facehi1,ilame_facehi2,ilame_facehi3
      REAL*8 lame_face(
     & ilame_facelo0:ilame_facehi0,
     & ilame_facelo1:ilame_facehi1,
     & ilame_facelo2:ilame_facehi2,
     & ilame_facelo3:ilame_facehi3,
     & 0:nlame_facecomp-1)
      integer nD_cfgcomp
      integer iD_cfglo0,iD_cfglo1,iD_cfglo2,iD_cfglo3
      integer iD_cfghi0,iD_cfghi1,iD_cfghi2,iD_cfghi3
      REAL*8 D_cfg(
     & iD_cfglo0:iD_cfghi0,
     & iD_cfglo1:iD_cfghi1,
     & iD_cfglo2:iD_cfghi2,
     & iD_cfglo3:iD_cfghi3,
     & 0:nD_cfgcomp-1)
      REAL*8 mass
      integer Nr
      integer simpleDiffusion
      integer ifBlo0,ifBlo1,ifBlo2,ifBlo3
      integer ifBhi0,ifBhi1,ifBhi2,ifBhi3
      REAL*8 fB(
     & ifBlo0:ifBhi0,
     & ifBlo1:ifBhi1,
     & ifBlo2:ifBhi2,
     & ifBlo3:ifBhi3)
      integer idfBdmulo0,idfBdmulo1,idfBdmulo2,idfBdmulo3
      integer idfBdmuhi0,idfBdmuhi1,idfBdmuhi2,idfBdmuhi3
      REAL*8 dfBdmu(
     & idfBdmulo0:idfBdmuhi0,
     & idfBdmulo1:idfBdmuhi1,
     & idfBdmulo2:idfBdmuhi2,
     & idfBdmulo3:idfBdmuhi3)
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL*8 B(
     & iBlo0:iBhi0,
     & iBlo1:iBhi1,
     & iBlo2:iBhi2,
     & iBlo3:iBhi3)
      integer iNlo0,iNlo1,iNlo2,iNlo3
      integer iNhi0,iNhi1,iNhi2,iNhi3
      REAL*8 N(
     & iNlo0:iNhi0,
     & iNlo1:iNhi1,
     & iNlo2:iNhi2,
     & iNlo3:iNhi3)
      integer iUlo0,iUlo1,iUlo2,iUlo3
      integer iUhi0,iUhi1,iUhi2,iUhi3
      REAL*8 U(
     & iUlo0:iUhi0,
     & iUlo1:iUhi1,
     & iUlo2:iUhi2,
     & iUlo3:iUhi3)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL*8 T(
     & iTlo0:iThi0,
     & iTlo1:iThi1,
     & iTlo2:iThi2,
     & iTlo3:iThi3)
      integer iClo0,iClo1,iClo2,iClo3
      integer iChi0,iChi1,iChi2,iChi3
      REAL*8 C(
     & iClo0:iChi0,
     & iClo1:iChi1,
     & iClo2:iChi2,
     & iClo3:iChi3)
      integer iPlo0,iPlo1,iPlo2,iPlo3
      integer iPhi0,iPhi1,iPhi2,iPhi3
      REAL*8 P(
     & iPlo0:iPhi0,
     & iPlo1:iPhi1,
     & iPlo2:iPhi2,
     & iPlo3:iPhi3)
      integer ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
      integer ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & ifluxlo3:ifluxhi3)
      integer i0,i1,i2,i3, L3, L4, kk
      double precision hr, htheta, hphi, vparr2, vperp2, v2, coef
      double precision fB_face, B_face, T_face, N_face, U_face, C_face
      double precision dfB_dr, dfB_dmu, dB_dr, dT_dr, dn_dr
      double precision FluxR, FluxZ, Upsi, Dpsi
      double precision D0, D1, D2, D3, DN
      double precision vp_RealCoord, mu_RealCoord
      kk = iBlo2
      L3 = iBlo3
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
         hr = lame_face(i0,i1,kk,L3,0)
         htheta = lame_face(i0,i1,kk,L3,1)
         hphi = lame_face(i0,i1,kk,L3,2)
         D0 = D_cfg(i0,i1,kk,L3,0)
         D1 = D_cfg(i0,i1,kk,L3,1)
         D2 = D_cfg(i0,i1,kk,L3,2)
         D3 = D_cfg(i0,i1,kk,L3,3)
         DN = D_cfg(i0,i1,kk,L3,4)
         if (dir==0) then
            dfB_dr = (fB(i0,i1,i2,i3)-fB(i0-1,i1,i2,i3))/dx(0)
            fB_face = (fB(i0,i1,i2,i3)+fB(i0-1,i1,i2,i3))/(2.0d0)
            dB_dr = (B(i0,i1,kk,L3)-B(i0-1,i1,kk,L3))/dx(0)
            B_face = (B(i0,i1,kk,L3)+B(i0-1,i1,kk,L3))/(2.0d0)
            N_face = (N(i0,i1,kk,L3)+N(i0-1,i1,kk,L3))/(2.0d0)
            T_face = (T(i0,i1,kk,L3)+T(i0-1,i1,kk,L3))/(2.0d0)
            U_face = (U(i0,i1,kk,L3)+U(i0-1,i1,kk,L3))/(2.0d0)
            C_face = (C(i0,i1,kk,L3)+C(i0-1,i1,kk,L3))/(2.0d0)
            dfB_dmu = (dfBdmu(i0,i1,i2,i3)+dfBdmu(i0-1,i1,i2,i3))/(2.0d0
     &)
            dN_dr = (N(i0,i1,kk,L3)-N(i0-1,i1,kk,L3))/dx(0)
            dT_dr = (T(i0,i1,kk,L3)-T(i0-1,i1,kk,L3))/dx(0)
            coef = (2.0d0)*C_face/(3.0d0)-(3.0d0)/(2.0d0)
            vp_RealCoord = (i2+(0.500d0))*dx(2)
            mu_RealCoord = (i3+(0.500d0))*dx(3)
            vparr2 = (vp_RealCoord-U_face)**2
            vperp2 = mu_RealCoord*B_face
            v2 = (mass*vparr2+vperp2)/(2.0d0)/T_face
            Upsi = (DN+D2/coef*(v2-(3.0d0)/(2.0d0)))*dN_dr/N_face/hr
     & + (D1+D3/coef*(v2-(3.0d0)/(2.0d0)))*dT_dr/T_face/hr
            Dpsi = D0
            if (simpleDiffusion.eq.1) then
               fluxR = hphi*htheta*( Dpsi*dfB_dr/hr )
            else
               fluxR = hphi*htheta*( Dpsi*(dfB_dr/hr - ( mu_RealCoord*df
     &B_dmu
     & + fB_face )*dB_dr/B_face/hr ) - Upsi*fB_face )
            endif
         else
            fluxR = 0.0
         endif
         flux(i0,i1,i2,i3) = fluxR
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVAL_ANOM_FLUX_NOT_ALIGNED_4D(
     & dir
     & ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & ,dx
     & ,NJinv
     & ,iNJinvlo0,iNJinvlo1,iNJinvlo2,iNJinvlo3
     & ,iNJinvhi0,iNJinvhi1,iNJinvhi2,iNJinvhi3
     & ,nNJinvcomp
     & ,bunit
     & ,ibunitlo0,ibunitlo1,ibunitlo2,ibunitlo3
     & ,ibunithi0,ibunithi1,ibunithi2,ibunithi3
     & ,nbunitcomp
     & ,D_cfg
     & ,iD_cfglo0,iD_cfglo1,iD_cfglo2,iD_cfglo3
     & ,iD_cfghi0,iD_cfghi1,iD_cfghi2,iD_cfghi3
     & ,nD_cfgcomp
     & ,mass
     & ,Nr
     & ,simpleDiffusion
     & ,fB
     & ,ifBlo0,ifBlo1,ifBlo2,ifBlo3
     & ,ifBhi0,ifBhi1,ifBhi2,ifBhi3
     & ,dfBdmu
     & ,idfBdmulo0,idfBdmulo1,idfBdmulo2,idfBdmulo3
     & ,idfBdmuhi0,idfBdmuhi1,idfBdmuhi2,idfBdmuhi3
     & ,B
     & ,iBlo0,iBlo1,iBlo2,iBlo3
     & ,iBhi0,iBhi1,iBhi2,iBhi3
     & ,N
     & ,iNlo0,iNlo1,iNlo2,iNlo3
     & ,iNhi0,iNhi1,iNhi2,iNhi3
     & ,U
     & ,iUlo0,iUlo1,iUlo2,iUlo3
     & ,iUhi0,iUhi1,iUhi2,iUhi3
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & ,C
     & ,iClo0,iClo1,iClo2,iClo3
     & ,iChi0,iChi1,iChi2,iChi3
     & ,P
     & ,iPlo0,iPlo1,iPlo2,iPlo3
     & ,iPhi0,iPhi1,iPhi2,iPhi3
     & ,flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
     & ,ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
     & ,nfluxcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL*8 dx(0:3)
      integer nNJinvcomp
      integer iNJinvlo0,iNJinvlo1,iNJinvlo2,iNJinvlo3
      integer iNJinvhi0,iNJinvhi1,iNJinvhi2,iNJinvhi3
      REAL*8 NJinv(
     & iNJinvlo0:iNJinvhi0,
     & iNJinvlo1:iNJinvhi1,
     & iNJinvlo2:iNJinvhi2,
     & iNJinvlo3:iNJinvhi3,
     & 0:nNJinvcomp-1)
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1,ibunitlo2,ibunitlo3
      integer ibunithi0,ibunithi1,ibunithi2,ibunithi3
      REAL*8 bunit(
     & ibunitlo0:ibunithi0,
     & ibunitlo1:ibunithi1,
     & ibunitlo2:ibunithi2,
     & ibunitlo3:ibunithi3,
     & 0:nbunitcomp-1)
      integer nD_cfgcomp
      integer iD_cfglo0,iD_cfglo1,iD_cfglo2,iD_cfglo3
      integer iD_cfghi0,iD_cfghi1,iD_cfghi2,iD_cfghi3
      REAL*8 D_cfg(
     & iD_cfglo0:iD_cfghi0,
     & iD_cfglo1:iD_cfghi1,
     & iD_cfglo2:iD_cfghi2,
     & iD_cfglo3:iD_cfghi3,
     & 0:nD_cfgcomp-1)
      REAL*8 mass
      integer Nr
      integer simpleDiffusion
      integer ifBlo0,ifBlo1,ifBlo2,ifBlo3
      integer ifBhi0,ifBhi1,ifBhi2,ifBhi3
      REAL*8 fB(
     & ifBlo0:ifBhi0,
     & ifBlo1:ifBhi1,
     & ifBlo2:ifBhi2,
     & ifBlo3:ifBhi3)
      integer idfBdmulo0,idfBdmulo1,idfBdmulo2,idfBdmulo3
      integer idfBdmuhi0,idfBdmuhi1,idfBdmuhi2,idfBdmuhi3
      REAL*8 dfBdmu(
     & idfBdmulo0:idfBdmuhi0,
     & idfBdmulo1:idfBdmuhi1,
     & idfBdmulo2:idfBdmuhi2,
     & idfBdmulo3:idfBdmuhi3)
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL*8 B(
     & iBlo0:iBhi0,
     & iBlo1:iBhi1,
     & iBlo2:iBhi2,
     & iBlo3:iBhi3)
      integer iNlo0,iNlo1,iNlo2,iNlo3
      integer iNhi0,iNhi1,iNhi2,iNhi3
      REAL*8 N(
     & iNlo0:iNhi0,
     & iNlo1:iNhi1,
     & iNlo2:iNhi2,
     & iNlo3:iNhi3)
      integer iUlo0,iUlo1,iUlo2,iUlo3
      integer iUhi0,iUhi1,iUhi2,iUhi3
      REAL*8 U(
     & iUlo0:iUhi0,
     & iUlo1:iUhi1,
     & iUlo2:iUhi2,
     & iUlo3:iUhi3)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL*8 T(
     & iTlo0:iThi0,
     & iTlo1:iThi1,
     & iTlo2:iThi2,
     & iTlo3:iThi3)
      integer iClo0,iClo1,iClo2,iClo3
      integer iChi0,iChi1,iChi2,iChi3
      REAL*8 C(
     & iClo0:iChi0,
     & iClo1:iChi1,
     & iClo2:iChi2,
     & iClo3:iChi3)
      integer iPlo0,iPlo1,iPlo2,iPlo3
      integer iPhi0,iPhi1,iPhi2,iPhi3
      REAL*8 P(
     & iPlo0:iPhi0,
     & iPlo1:iPhi1,
     & iPlo2:iPhi2,
     & iPlo3:iPhi3)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
      integer ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & ifluxlo3:ifluxhi3,
     & 0:nfluxcomp-1)
      integer i0,i1,i2,i3, L3, L4, kk
      double precision fB_face, dfB_dpsi, dfB_dtheta, dfB_dmu
      double precision dB_dpsi, B_face, dB_dtheta
      double precision fluxR, fluxZ, dot_product, bpol
      double precision N_face, T_face, U_face, C_face
      double precision dN_dpsi, dT_dpsi, dN_dtheta, dT_dtheta
      double precision vparr2, vperp2, v2, coef, Upsi, Utheta
      double precision fluxpsi, fluxtheta
      double precision D0, D1, D2, D3, DN
      double precision vp_RealCoord, mu_RealCoord
      kk = iBlo2
      L3 = iBlo3
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
         D0 = D_cfg(i0,i1,kk,L3,0)
         D1 = D_cfg(i0,i1,kk,L3,1)
         D2 = D_cfg(i0,i1,kk,L3,2)
         D3 = D_cfg(i0,i1,kk,L3,3)
         DN = D_cfg(i0,i1,kk,L3,4)
         if (dir.eq.0) then
            fB_face = (fB(i0,i1,i2,i3)+fB(i0-1,i1,i2,i3))/(2.0d0)
            B_face = (B(i0,i1,kk,L3)+B(i0-1,i1,kk,L3))/(2.0d0)
            N_face = (N(i0,i1,kk,L3)+N(i0-1,i1,kk,L3))/(2.0d0)
            T_face = (T(i0,i1,kk,L3)+T(i0-1,i1,kk,L3))/(2.0d0)
            U_face = (U(i0,i1,kk,L3)+U(i0-1,i1,kk,L3))/(2.0d0)
            C_face = (C(i0,i1,kk,L3)+C(i0-1,i1,kk,L3))/(2.0d0)
            dfB_dpsi= (fB(i0,i1,i2,i3)-fB(i0-1,i1,i2,i3))/dx(0)
            dB_dpsi = (B(i0,i1,kk,L3)-B(i0-1,i1,kk,L3))/dx(0)
            dN_dpsi = (N(i0,i1,kk,L3)-N(i0-1,i1,kk,L3))/dx(0)
            dT_dpsi = (T(i0,i1,kk,L3)-T(i0-1,i1,kk,L3))/dx(0)
            dfB_dtheta= (fB(i0,i1+1,i2,i3) +fB(i0-1,i1+1,i2,i3) -fB(i0,i
     &1-1,i2,i3) -fB(i0-1,i1-1,i2,i3) )/(4.0d0)/dx(1)
            dB_dtheta = (B(i0,i1+1,kk,L3)+B(i0-1,i1+1,kk,L3)-B(i0,i1-1,k
     &k,L3)-B(i0-1,i1-1,kk,L3))/(4.0d0)/dx(1)
            dN_dtheta = (N(i0,i1+1,kk,L3)+N(i0-1,i1+1,kk,L3)-N(i0,i1-1,k
     &k,L3)-N(i0-1,i1-1,kk,L3))/(4.0d0)/dx(1)
            dT_dtheta = (T(i0,i1+1,kk,L3)+T(i0-1,i1+1,kk,L3)-T(i0,i1-1,k
     &k,L3)-T(i0-1,i1-1,kk,L3))/(4.0d0)/dx(1)
            dfB_dmu = (dfBdmu(i0,i1,i2,i3)+dfBdmu(i0-1,i1,i2,i3))/(2.0d0
     &)
         else if (dir.eq.1) then
            fB_face = (fB(i0,i1,i2,i3)+fB(i0,i1-1,i2,i3))/(2.0d0)
            B_face = (B(i0,i1,kk,L3)+B(i0,i1-1,kk,L3))/(2.0d0)
            N_face = (N(i0,i1,kk,L3)+N(i0,i1-1,kk,L3))/(2.0d0)
            T_face = (T(i0,i1,kk,L3)+T(i0,i1-1,kk,L3))/(2.0d0)
            U_face = (U(i0,i1,kk,L3)+U(i0,i1-1,kk,L3))/(2.0d0)
            C_face = (C(i0,i1,kk,L3)+C(i0,i1-1,kk,L3))/(2.0d0)
            dfB_dtheta= (fB(i0,i1,i2,i3)-fB(i0,i1-1,i2,i3))/dx(1)
            dB_dtheta = (B(i0,i1,kk,L3)-B(i0,i1-1,kk,L3))/dx(1)
            dN_dtheta = (N(i0,i1,kk,L3)-N(i0,i1-1,kk,L3))/dx(1)
            dT_dtheta = (T(i0,i1,kk,L3)-T(i0,i1-1,kk,L3))/dx(1)
            dfB_dpsi= (fB(i0+1,i1,i2,i3) +fB(i0+1,i1-1,i2,i3) -fB(i0-1,i
     &1,i2,i3) -fB(i0-1,i1-1,i2,i3) )/(4.0d0)/dx(0)
            dB_dpsi = (B(i0+1,i1,kk,L3)+B(i0+1,i1-1,kk,L3)-B(i0-1,i1,kk,
     &L3)-B(i0-1,i1-1,kk,L3))/(4.0d0)/dx(0)
            dN_dpsi = (N(i0+1,i1,kk,L3)+N(i0+1,i1-1,kk,L3)-N(i0-1,i1,kk,
     &L3)-N(i0-1,i1-1,kk,L3))/(4.0d0)/dx(0)
            dT_dpsi = (T(i0+1,i1,kk,L3)+T(i0+1,i1-1,kk,L3)-T(i0-1,i1,kk,
     &L3)-T(i0-1,i1-1,kk,L3))/(4.0d0)/dx(0)
            dfB_dmu = (dfBdmu(i0,i1,i2,i3)+dfBdmu(i0,i1-1,i2,i3))/(2.0d0
     &)
         else
            dfB_dpsi = (fB(i0+1,i1,i2,i3)-fB(i0-1,i1,i2,i3))/(2.0d0)/dx(
     &0)
            dfB_dtheta = (fB(i0,i1+1,i2,i3)-fB(i0,i1-1,i2,i3))/(2.0d0)/d
     &x(1)
         endif
         coef = (2.0d0)*C_face/(3.0d0)-(3.0d0)/(2.0d0)
         vp_RealCoord = (i2+(0.500d0))*dx(2)
         mu_RealCoord = (i3+(0.500d0))*dx(3)
         vparr2 = (vp_RealCoord-U_face)**2
         vperp2 = mu_RealCoord*B_face
         v2 = (mass*vparr2+vperp2)/(2.0d0)/T_face
         Upsi = (DN+D2/coef*(v2-(3.0d0)/(2.0d0)))*dN_dpsi/N_face
     & + (D1+D3/coef*(v2-(3.0d0)/(2.0d0)))*dT_dpsi/T_face
         Utheta = (DN+D2/coef*(v2-(3.0d0)/(2.0d0)))*dN_dtheta/N_face
     & + (D1+D3/coef*(v2-(3.0d0)/(2.0d0)))*dT_dtheta/T_face
         if (simpleDiffusion.eq.1) then
            fluxpsi = D0*dfB_dpsi
            fluxtheta = D0*dfB_dtheta
         else
            fluxpsi = D0*(dfB_dpsi - ( mu_RealCoord*dfB_dmu + fB_face )*
     &dB_dpsi/B_face)
     & - Upsi*fB_face
            fluxtheta = D0*(dfB_dtheta - ( mu_RealCoord*dfB_dmu + fB_fac
     &e )*dB_dtheta/B_face)
     & - Utheta*fB_face
         endif
         fluxR = NJinv(i0,i1,kk,L3,0) * fluxpsi + NJinv(i0,i1,kk,L3,1) *
     & fluxtheta
         fluxZ = NJinv(i0,i1,kk,L3,2) * fluxpsi + NJinv(i0,i1,kk,L3,3) *
     & fluxtheta
         bpol = sqrt(bunit(i0,i1,kk,L3,0)**2 + bunit(i0,i1,kk,L3,2)**2)
         dot_product = fluxR * bunit(i0,i1,kk,L3,2)/bpol - fluxZ * bunit
     &(i0,i1,kk,L3,0)/bpol
         fluxR = bunit(i0,i1,kk,L3,2)/bpol * dot_product
         fluxZ = -bunit(i0,i1,kk,L3,0)/bpol * dot_product
         flux(i0,i1,i2,i3,0) = fluxR
         flux(i0,i1,i2,i3,1) = fluxZ
         flux(i0,i1,i2,i3,2) = (0.0d0)
         flux(i0,i1,i2,i3,3) = (0.0d0)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine DFBDMU_CELL_CENTER_4D(
     & igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & ,dx
     & ,nmu
     & ,fB
     & ,ifBlo0,ifBlo1,ifBlo2,ifBlo3
     & ,ifBhi0,ifBhi1,ifBhi2,ifBhi3
     & ,dfB_dmu
     & ,idfB_dmulo0,idfB_dmulo1,idfB_dmulo2,idfB_dmulo3
     & ,idfB_dmuhi0,idfB_dmuhi1,idfB_dmuhi2,idfB_dmuhi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL*8 dx(0:3)
      integer nmu
      integer ifBlo0,ifBlo1,ifBlo2,ifBlo3
      integer ifBhi0,ifBhi1,ifBhi2,ifBhi3
      REAL*8 fB(
     & ifBlo0:ifBhi0,
     & ifBlo1:ifBhi1,
     & ifBlo2:ifBhi2,
     & ifBlo3:ifBhi3)
      integer idfB_dmulo0,idfB_dmulo1,idfB_dmulo2,idfB_dmulo3
      integer idfB_dmuhi0,idfB_dmuhi1,idfB_dmuhi2,idfB_dmuhi3
      REAL*8 dfB_dmu(
     & idfB_dmulo0:idfB_dmuhi0,
     & idfB_dmulo1:idfB_dmuhi1,
     & idfB_dmulo2:idfB_dmuhi2,
     & idfB_dmulo3:idfB_dmuhi3)
      integer i,j,k,l
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
       if (l==0) then
         dfB_dmu(i,j,k,l) = ( (0.0d0) + ( fB(i,j,k,l+1)-fB(i,j,k,l) )/dx
     &(3) )/(2.0d0)
       else
       if (l==nmu-1) then
         dfB_dmu(i,j,k,l) = ( (3.0d0)*fB(i,j,k,l)-(4.0d0)*fB(i,j,k,l-1)+
     &fB(i,j,k,l-2) )/(2.0d0)/dx(3)
       else
         dfB_dmu(i,j,k,l) = ( fB(i,j,k,l+1)-fB(i,j,k,l-1) )/(2.0d0)/dx(3
     &)
       endif
       endif
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine LAME_COEFFICIENTS_4D(
     & igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & ,N
     & ,iNlo0,iNlo1,iNlo2,iNlo3
     & ,iNhi0,iNhi1,iNhi2,iNhi3
     & ,nNcomp
     & ,dXdq
     & ,idXdqlo0,idXdqlo1,idXdqlo2,idXdqlo3
     & ,idXdqhi0,idXdqhi1,idXdqhi2,idXdqhi3
     & ,ndXdqcomp
     & ,Lame_coeff
     & ,iLame_coefflo0,iLame_coefflo1,iLame_coefflo2,iLame_coefflo3
     & ,iLame_coeffhi0,iLame_coeffhi1,iLame_coeffhi2,iLame_coeffhi3
     & ,nLame_coeffcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer nNcomp
      integer iNlo0,iNlo1,iNlo2,iNlo3
      integer iNhi0,iNhi1,iNhi2,iNhi3
      REAL*8 N(
     & iNlo0:iNhi0,
     & iNlo1:iNhi1,
     & iNlo2:iNhi2,
     & iNlo3:iNhi3,
     & 0:nNcomp-1)
      integer ndXdqcomp
      integer idXdqlo0,idXdqlo1,idXdqlo2,idXdqlo3
      integer idXdqhi0,idXdqhi1,idXdqhi2,idXdqhi3
      REAL*8 dXdq(
     & idXdqlo0:idXdqhi0,
     & idXdqlo1:idXdqhi1,
     & idXdqlo2:idXdqhi2,
     & idXdqlo3:idXdqhi3,
     & 0:ndXdqcomp-1)
      integer nLame_coeffcomp
      integer iLame_coefflo0,iLame_coefflo1,iLame_coefflo2,iLame_coefflo
     &3
      integer iLame_coeffhi0,iLame_coeffhi1,iLame_coeffhi2,iLame_coeffhi
     &3
      REAL*8 Lame_coeff(
     & iLame_coefflo0:iLame_coeffhi0,
     & iLame_coefflo1:iLame_coeffhi1,
     & iLame_coefflo2:iLame_coeffhi2,
     & iLame_coefflo3:iLame_coeffhi3,
     & 0:nLame_coeffcomp-1)
      integer i0,i1,i2,i3, L3, L4, kk
      double precision hr, htheta, hphi, Jacobian
      kk = iNlo2
      L3 = iNlo3
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
        hr = sqrt(dXdq(i0,i1,kk,L3,0)**2 + dXdq(i0,i1,kk,L3,2)**2)
        htheta = sqrt(dXdq(i0,i1,kk,L3,1)**2 + dXdq(i0,i1,kk,L3,3)**2)
        Jacobian = N(i0,i1,kk,L3,0)*dXdq(i0,i1,kk,L3,0) + N(i0,i1,kk,L3,
     &1)*dXdq(i0,i1,kk,L3,1)
        hphi = Jacobian/hr/htheta
        Lame_coeff(i0,i1,i2,i3,0) = hr
        Lame_coeff(i0,i1,i2,i3,1) = htheta
        Lame_coeff(i0,i1,i2,i3,2) = hphi
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EVAL_BETA_4D(
     & igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     & ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     & ,dx
     & ,D_kinet
     & ,iD_kinethi0
     & ,mass
     & ,Nr
     & ,B
     & ,iBlo0,iBlo1,iBlo2,iBlo3
     & ,iBhi0,iBhi1,iBhi2,iBhi3
     & ,N
     & ,iNlo0,iNlo1,iNlo2,iNlo3
     & ,iNhi0,iNhi1,iNhi2,iNhi3
     & ,U
     & ,iUlo0,iUlo1,iUlo2,iUlo3
     & ,iUhi0,iUhi1,iUhi2,iUhi3
     & ,T
     & ,iTlo0,iTlo1,iTlo2,iTlo3
     & ,iThi0,iThi1,iThi2,iThi3
     & ,C
     & ,iClo0,iClo1,iClo2,iClo3
     & ,iChi0,iChi1,iChi2,iChi3
     & ,P
     & ,iPlo0,iPlo1,iPlo2,iPlo3
     & ,iPhi0,iPhi1,iPhi2,iPhi3
     & ,beta
     & ,ibetalo0,ibetalo1,ibetalo2,ibetalo3
     & ,ibetahi0,ibetahi1,ibetahi2,ibetahi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL*8 dx(0:3)
      integer iD_kinethi0
      REAL*8 D_kinet(
     & 0:iD_kinethi0)
      REAL*8 mass
      integer Nr
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL*8 B(
     & iBlo0:iBhi0,
     & iBlo1:iBhi1,
     & iBlo2:iBhi2,
     & iBlo3:iBhi3)
      integer iNlo0,iNlo1,iNlo2,iNlo3
      integer iNhi0,iNhi1,iNhi2,iNhi3
      REAL*8 N(
     & iNlo0:iNhi0,
     & iNlo1:iNhi1,
     & iNlo2:iNhi2,
     & iNlo3:iNhi3)
      integer iUlo0,iUlo1,iUlo2,iUlo3
      integer iUhi0,iUhi1,iUhi2,iUhi3
      REAL*8 U(
     & iUlo0:iUhi0,
     & iUlo1:iUhi1,
     & iUlo2:iUhi2,
     & iUlo3:iUhi3)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL*8 T(
     & iTlo0:iThi0,
     & iTlo1:iThi1,
     & iTlo2:iThi2,
     & iTlo3:iThi3)
      integer iClo0,iClo1,iClo2,iClo3
      integer iChi0,iChi1,iChi2,iChi3
      REAL*8 C(
     & iClo0:iChi0,
     & iClo1:iChi1,
     & iClo2:iChi2,
     & iClo3:iChi3)
      integer iPlo0,iPlo1,iPlo2,iPlo3
      integer iPhi0,iPhi1,iPhi2,iPhi3
      REAL*8 P(
     & iPlo0:iPhi0,
     & iPlo1:iPhi1,
     & iPlo2:iPhi2,
     & iPlo3:iPhi3)
      integer ibetalo0,ibetalo1,ibetalo2,ibetalo3
      integer ibetahi0,ibetahi1,ibetahi2,ibetahi3
      REAL*8 beta(
     & ibetalo0:ibetahi0,
     & ibetalo1:ibetahi1,
     & ibetalo2:ibetahi2,
     & ibetalo3:ibetahi3)
      integer i0,i1,i2,i3, L3, L4, kk
      double precision vparr2, vperp2, v2, coef, dlnB_dr, dlnN_dr, dlnT_
     &dr
      double precision DB, Upsi
      double precision vp_RealCoord, mu_RealCoord
      kk = iBlo2
      L3 = iBlo3
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
       DB = -(2.0d0)/(3.0d0)*D_kinet(0)*P(i0,i1,kk,L3)
       dlnB_dr = (B(i0+1,i1,kk,L3)-B(i0-1,i1,kk,L3))/(2.0d0)/dx(0)/B(i0,
     &i1,kk,L3)
       if (i0==0) then
         dlnN_dr = (-(3.0d0)*N(i0,i1,kk,L3)+(4.0d0)*N(i0+1,i1,kk,L3)-N(i
     &0+2,i1,kk,L3))/(2.0d0)/dx(0)/N(i0,i1,kk,L3)
         dlnT_dr = (-(3.0d0)*T(i0,i1,kk,L3)+(4.0d0)*T(i0+1,i1,kk,L3)-T(i
     &0+2,i1,kk,L3))/(2.0d0)/dx(0)/T(i0,i1,kk,L3)
       else
       if (i0==Nr-1) then
         dlnN_dr = ((3.0d0)*N(i0,i1,kk,L3)-(4.0d0)*N(i0-1,i1,kk,L3)+N(i0
     &-2,i1,kk,L3))/(2.0d0)/dx(0)/N(i0,i1,kk,L3)
         dlnT_dr = ((3.0d0)*T(i0,i1,kk,L3)-(4.0d0)*T(i0-1,i1,kk,L3)+T(i0
     &-2,i1,kk,L3))/(2.0d0)/dx(0)/T(i0,i1,kk,L3)
       else
         dlnN_dr = (N(i0+1,i1,kk,L3)-N(i0-1,i1,kk,L3))/(2.0d0)/dx(0)/N(i
     &0,i1,kk,L3)
         dlnT_dr = (T(i0+1,i1,kk,L3)-T(i0-1,i1,kk,L3))/(2.0d0)/dx(0)/T(i
     &0,i1,kk,L3)
       endif
       endif
       coef = (2.0d0)/(3.0d0)*C(i0,i1,kk,L3)-(3.0d0)/(2.0d0)
       vp_RealCoord = (i2+(0.500d0))*dx(2)
       mu_RealCoord = (i3+(0.500d0))*dx(3)
       vparr2 = (vp_RealCoord-U(i0,i1,kk,L3))**2
       vperp2 = mu_RealCoord*B(i0,i1,kk,L3)
       v2 = (mass*vparr2+vperp2)/(2.0d0)/T(i0,i1,kk,L3)
       Upsi = D_kinet(2)/coef*(v2-(3.0d0)/(2.0d0))*dlnN_dr + (D_kinet(1)
     &+D_kinet(3)/coef*(v2-(3.0d0)/(2.0d0)))*dlnT_dr
       Upsi = Upsi + DB*coef*(v2-(3.0d0)/(2.0d0))*dlnB_dr
       beta(i0,i1,i2,i3) = Upsi**2/D_kinet(0)
      enddo
      enddo
      enddo
      enddo
      return
      end
