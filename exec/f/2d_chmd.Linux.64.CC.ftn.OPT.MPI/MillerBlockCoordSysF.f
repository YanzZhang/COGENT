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
      subroutine GET_MILLER_FIELD_DATA_2D(
     & igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,RBpol
     & ,iRBpollo0,iRBpollo1
     & ,iRBpolhi0,iRBpolhi1
     & ,dRBpoldt
     & ,idRBpoldtlo0,idRBpoldtlo1
     & ,idRBpoldthi0,idRBpoldthi1
     & ,RBtor
     & ,Rmaj
     & ,iRmajlo0,iRmajlo1
     & ,iRmajhi0,iRmajhi1
     & ,Rr
     & ,iRrlo0,iRrlo1
     & ,iRrhi0,iRrhi1
     & ,Rt
     & ,iRtlo0,iRtlo1
     & ,iRthi0,iRthi1
     & ,Zr
     & ,iZrlo0,iZrlo1
     & ,iZrhi0,iZrhi1
     & ,Zt
     & ,iZtlo0,iZtlo1
     & ,iZthi0,iZthi1
     & ,Rtt
     & ,iRttlo0,iRttlo1
     & ,iRtthi0,iRtthi1
     & ,Ztt
     & ,iZttlo0,iZttlo1
     & ,iZtthi0,iZtthi1
     & ,b_pt
     & ,ib_ptlo0,ib_ptlo1
     & ,ib_pthi0,ib_pthi1
     & ,nb_ptcomp
     & ,Bmag_pt
     & ,iBmag_ptlo0,iBmag_ptlo1
     & ,iBmag_pthi0,iBmag_pthi1
     & ,bunit_pt
     & ,ibunit_ptlo0,ibunit_ptlo1
     & ,ibunit_pthi0,ibunit_pthi1
     & ,nbunit_ptcomp
     & ,gradb_pt
     & ,igradb_ptlo0,igradb_ptlo1
     & ,igradb_pthi0,igradb_pthi1
     & ,ngradb_ptcomp
     & ,curlb_pt
     & ,icurlb_ptlo0,icurlb_ptlo1
     & ,icurlb_pthi0,icurlb_pthi1
     & ,ncurlb_ptcomp
     & ,bdotcurlb_pt
     & ,ibdotcurlb_ptlo0,ibdotcurlb_ptlo1
     & ,ibdotcurlb_pthi0,ibdotcurlb_pthi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer iRBpollo0,iRBpollo1
      integer iRBpolhi0,iRBpolhi1
      REAL*8 RBpol(
     & iRBpollo0:iRBpolhi0,
     & iRBpollo1:iRBpolhi1)
      integer idRBpoldtlo0,idRBpoldtlo1
      integer idRBpoldthi0,idRBpoldthi1
      REAL*8 dRBpoldt(
     & idRBpoldtlo0:idRBpoldthi0,
     & idRBpoldtlo1:idRBpoldthi1)
      REAL*8 RBtor
      integer iRmajlo0,iRmajlo1
      integer iRmajhi0,iRmajhi1
      REAL*8 Rmaj(
     & iRmajlo0:iRmajhi0,
     & iRmajlo1:iRmajhi1)
      integer iRrlo0,iRrlo1
      integer iRrhi0,iRrhi1
      REAL*8 Rr(
     & iRrlo0:iRrhi0,
     & iRrlo1:iRrhi1)
      integer iRtlo0,iRtlo1
      integer iRthi0,iRthi1
      REAL*8 Rt(
     & iRtlo0:iRthi0,
     & iRtlo1:iRthi1)
      integer iZrlo0,iZrlo1
      integer iZrhi0,iZrhi1
      REAL*8 Zr(
     & iZrlo0:iZrhi0,
     & iZrlo1:iZrhi1)
      integer iZtlo0,iZtlo1
      integer iZthi0,iZthi1
      REAL*8 Zt(
     & iZtlo0:iZthi0,
     & iZtlo1:iZthi1)
      integer iRttlo0,iRttlo1
      integer iRtthi0,iRtthi1
      REAL*8 Rtt(
     & iRttlo0:iRtthi0,
     & iRttlo1:iRtthi1)
      integer iZttlo0,iZttlo1
      integer iZtthi0,iZtthi1
      REAL*8 Ztt(
     & iZttlo0:iZtthi0,
     & iZttlo1:iZtthi1)
      integer nb_ptcomp
      integer ib_ptlo0,ib_ptlo1
      integer ib_pthi0,ib_pthi1
      REAL*8 b_pt(
     & ib_ptlo0:ib_pthi0,
     & ib_ptlo1:ib_pthi1,
     & 0:nb_ptcomp-1)
      integer iBmag_ptlo0,iBmag_ptlo1
      integer iBmag_pthi0,iBmag_pthi1
      REAL*8 Bmag_pt(
     & iBmag_ptlo0:iBmag_pthi0,
     & iBmag_ptlo1:iBmag_pthi1)
      integer nbunit_ptcomp
      integer ibunit_ptlo0,ibunit_ptlo1
      integer ibunit_pthi0,ibunit_pthi1
      REAL*8 bunit_pt(
     & ibunit_ptlo0:ibunit_pthi0,
     & ibunit_ptlo1:ibunit_pthi1,
     & 0:nbunit_ptcomp-1)
      integer ngradb_ptcomp
      integer igradb_ptlo0,igradb_ptlo1
      integer igradb_pthi0,igradb_pthi1
      REAL*8 gradb_pt(
     & igradb_ptlo0:igradb_pthi0,
     & igradb_ptlo1:igradb_pthi1,
     & 0:ngradb_ptcomp-1)
      integer ncurlb_ptcomp
      integer icurlb_ptlo0,icurlb_ptlo1
      integer icurlb_pthi0,icurlb_pthi1
      REAL*8 curlb_pt(
     & icurlb_ptlo0:icurlb_pthi0,
     & icurlb_ptlo1:icurlb_pthi1,
     & 0:ncurlb_ptcomp-1)
      integer ibdotcurlb_ptlo0,ibdotcurlb_ptlo1
      integer ibdotcurlb_pthi0,ibdotcurlb_pthi1
      REAL*8 bdotcurlb_pt(
     & ibdotcurlb_ptlo0:ibdotcurlb_pthi0,
     & ibdotcurlb_ptlo1:ibdotcurlb_pthi1)
      integer i,j, l
      double precision
     & jac, RB, bdotcurlb, bunit(0:2),
     & bunitR, dbunitRdt, bunitphi, dbunitphidt, bunitZ, dbunitZdt,
     & dRBdt, rsphert, rsphertsq, wronsk, curlb(0:2), bporb
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         RB = dsqrt(RBpol(i,j)**2 + RBtor**2)
         jac = Rr(i,j) * Zt(i,j) - Zr(i,j) * Rt(i,j)
         dRBdt = RBpol(i,j) * dRBpoldt(i,j) / RB
         rsphertsq = Zt(i,j)**2 + Rt(i,j)**2
         rsphert = dsqrt(rsphertsq)
         bunitR = RBpol(i,j) * Rt(i,j) / (rsphert*RB)
         bunitphi = RBtor / RB
         bunitZ = RBpol(i,j) * Zt(i,j) / (rsphert*RB)
         bunit(0) = bunitR
         bunit(1) = bunitphi
         bunit(2) = bunitZ
         wronsk = (Rt(i,j)*Ztt(i,j) - Zt(i,j)*Rtt(i,j))/rsphertsq
         dbunitRdt = ((dRBpoldt(i,j)/RBpol(i,j))-drBdt/RB)*bunitR - buni
     &tZ*wronsk
         dbunitZdt = ((dRBpoldt(i,j)/RBpol(i,j))-drBdt/RB)*bunitZ + buni
     &tZ*wronsk
         dbunitphidt = -(RBtor/(RB)**2)*dRBdt
         do l = 0, 2
            bunit_pt(i,j,l) = bunit(l)
            b_pt(i,j,l) = bunit(l) * (RB/Rmaj(i,j))
         end do
         bporb = RBpol(i,j)/(Rmaj(i,j)*RB*jac)
         gradb_pt(i,j,0) = - (RB/Rmaj(i,j)**2)-Zr(i,j)*bporb*dRBdt
         gradb_pt(i,j,1) = (0.0d0)
         gradb_pt(i,j,2) = Rr(i,j)*bporb*dRBdt
         curlb(0) = -Rr(i,j)*dbunitphidt/jac
         curlb(1) = (Rr(i,j)*dbunitRdt + Zr(i,j)*dbunitZdt)/jac
         curlb(2) = bunitphi/Rmaj(i,j) - Zr(i,j)*dbunitphidt/jac
         bdotcurlb = (0.0d0)
         do l=0, 2
            curlb_pt(i,j,l) = curlb(l)
            bdotcurlb = bdotcurlb + curlb(l) * bunit(l)
         end do
         Bmag_pt(i,j) = RB/Rmaj(i,j)
         bdotcurlb_pt(i,j) = bdotcurlb
      enddo
      enddo
      return
      end
      function BpolR(t, beta, kappa, dpsidr, drr0, sk, sd)
      implicit none
      double precision BpolR, t, beta, kappa, dpsidr, drr0, sk, sd
      double precision st, ct, tpbs
      st = dsin(t)
      ct = dcos(t)
      tpbs = t + beta*st
      BpolR = (dpsidr / kappa) *
     & dsqrt( (dsin(tpbs)*((1.0d0) + beta*ct))**2 + (kappa*ct)**2) /
     & (dcos(beta*st) + drr0*ct + (sk - sd*ct + ((1.0d0) + sk)*beta*ct)
     & *st*dsin(tpbs))
      return
      end
      function dBpolRdt(t, beta, kappa, dpsidr, drr0, sk, sd)
      implicit none
      double precision dBpolRdt, t, beta, kappa, dpsidr, drr0, sk, sd
      double precision st, ct, tpbs, stpbs, ctpbs, opct, n, d, dndt, ddd
     &t
      st = dsin(t)
      ct = dcos(t)
      tpbs = t + beta*st
      stpbs = dsin(tpbs)
      ctpbs = dcos(tpbs)
      opct = (1.0d0) + beta*ct
      n = dsqrt( (stpbs*opct)**2 + (kappa*ct)**2)
      d = dcos(beta*st) + drr0*ct + (sk - sd*ct
     & + ((1.0d0) + sk)*beta*ct)*st*dsin(tpbs)
      dndt = ( opct*stpbs*(-beta*stpbs*st + opct**2*ctpbs ) - kappa**2 *
     & st*ct )/n
      dddt = -dsin(beta*st)*beta*ct - drr0*st
     & + (sk - sd*ct + ((1.0d0) + sk)*beta*ct)*(st*ctpbs*opct + ct*stpbs
     &)
     & + (sd - ((1.0d0) + sk)*beta) * st*st * stpbs
      dBpolRdt = (dpsidr / kappa) * (d*dndt - n*dddt) / d**2
      return
      end
      subroutine GET_RBPOL_MILLER_2D(
     & dir
     & ,igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,h
     & ,beta
     & ,kappa
     & ,dpsidr
     & ,r0
     & ,drr0
     & ,s_kappa
     & ,s_delta
     & ,RBpol
     & ,iRBpollo0,iRBpollo1
     & ,iRBpolhi0,iRBpolhi1
     & ,dRBpoldt
     & ,idRBpoldtlo0,idRBpoldtlo1
     & ,idRBpoldthi0,idRBpoldthi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer dir
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      REAL*8 h(0:1)
      REAL*8 beta
      REAL*8 kappa
      REAL*8 dpsidr
      REAL*8 r0
      REAL*8 drr0
      REAL*8 s_kappa
      REAL*8 s_delta
      integer iRBpollo0,iRBpollo1
      integer iRBpolhi0,iRBpolhi1
      REAL*8 RBpol(
     & iRBpollo0:iRBpolhi0,
     & iRBpollo1:iRBpolhi1)
      integer idRBpoldtlo0,idRBpoldtlo1
      integer idRBpoldthi0,idRBpoldthi1
      REAL*8 dRBpoldt(
     & idRBpoldtlo0:idRBpoldthi0,
     & idRBpoldtlo1:idRBpoldthi1)
      integer i,j, l
      double precision t, BpolR, dBpolRdt
      external BpolR, dBpolRdt
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         t = ( j + (0.500d0)*(1-CHF_ID(1,dir)) )*h(1)
         RBpol(i,j) = BpolR(t, beta, kappa, dpsidr, drr0, s_kappa, s_del
     &ta)
         dRBpoldt(i,j) = dBpolRdt(t, beta, kappa, dpsidr, drr0, s_kappa,
     & s_delta)
      enddo
      enddo
      return
      end
      subroutine MILLER_DXDXI_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,beta
     & ,kappa
     & ,Xi
     & ,iXilo0,iXilo1
     & ,iXihi0,iXihi1
     & ,nXicomp
     & ,destComp
     & ,dirX
     & ,dirXi
     & ,dXdXi
     & ,idXdXilo0,idXdXilo1
     & ,idXdXihi0,idXdXihi1
     & ,ndXdXicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 beta
      REAL*8 kappa
      integer nXicomp
      integer iXilo0,iXilo1
      integer iXihi0,iXihi1
      REAL*8 Xi(
     & iXilo0:iXihi0,
     & iXilo1:iXihi1,
     & 0:nXicomp-1)
      integer destComp
      integer dirX
      integer dirXi
      integer ndXdXicomp
      integer idXdXilo0,idXdXilo1
      integer idXdXihi0,idXdXihi1
      REAL*8 dXdXi(
     & idXdXilo0:idXdXihi0,
     & idXdXilo1:idXdXihi1,
     & 0:ndXdXicomp-1)
      integer i,j
      double precision r, theta
      if (dirX .eq. 0) then
         if (dirXi .eq. 0) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
               r = Xi(i,j,0)
               theta = Xi(i,j,1)
               dXdXi(i,j,destComp) = dcos(theta + beta*dsin(theta))
      enddo
      enddo
         else if (dirXi == 1) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
               r = Xi(i,j,0)
               theta = Xi(i,j,1)
               dXdXi(i,j,destComp) = -r*dsin(theta + beta*dsin(theta))*(
     &1. + beta*dcos(theta))
      enddo
      enddo
         endif
      else if (dirX .eq. 1) then
         if (dirXi .eq. 0) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
               r = Xi(i,j,0)
               theta = Xi(i,j,1)
               dXdXi(i,j,destComp) = kappa * dsin(theta)
      enddo
      enddo
         else if (dirXi .eq. 1) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
               r = Xi(i,j,0)
               theta = Xi(i,j,1)
               dXdXi(i,j,destComp) = r * kappa * dcos(theta)
      enddo
      enddo
         endif
      endif
      return
      end
      subroutine GET_MILLER_NC_MAPPED_COORDS_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dx
     & ,rmid
     & ,constminorrad
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 dx(0:1)
      REAL*8 rmid
      integer constminorrad
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      if (constminorrad .ne. 0) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            xi(i,j,0) = rmid;
            xi(i,j,1) = j*dx(1);
      enddo
      enddo
      else
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            xi(i,j,0) = i*dx(0);
            xi(i,j,1) = j*dx(1);
      enddo
      enddo
      endif
      return
      end
      subroutine GET_MILLER_CC_MAPPED_COORDS_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dx
     & ,rmid
     & ,constminorrad
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 dx(0:1)
      REAL*8 rmid
      integer constminorrad
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      if (constminorrad .ne. 0) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            xi(i,j,0) = rmid;
            xi(i,j,1) = (j + (0.500d0))*dx(1);
      enddo
      enddo
      else
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            xi(i,j,0) = (i + (0.500d0))*dx(0);
            xi(i,j,1) = (j + (0.500d0))*dx(1);
      enddo
      enddo
      endif
      return
      end
      subroutine GET_MILLER_FC_MAPPED_COORDS_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dir
     & ,dx
     & ,rmid
     & ,constminorrad
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL*8 dx(0:1)
      REAL*8 rmid
      integer constminorrad
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      double precision offset(0:1)
      offset(0) = (0.500d0)
      offset(1) = (0.500d0)
      offset(dir) = (0.0d0)
      if (constminorrad .ne. 0) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            xi(i,j,0) = rmid;
            xi(i,j,1) = (j + offset(1))*dx(1);
      enddo
      enddo
      else
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            xi(i,j,0) = (i + offset(0))*dx(0);
            xi(i,j,1) = (j + offset(1))*dx(1);
      enddo
      enddo
      endif
      return
      end
