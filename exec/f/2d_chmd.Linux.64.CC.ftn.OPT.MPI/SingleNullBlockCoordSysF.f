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
      subroutine GET_FIELD_FROM_RBPOL_RBTOR_2D(
     & igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,RBtor
     & ,RBpol
     & ,iRBpollo0,iRBpollo1
     & ,iRBpolhi0,iRBpolhi1
     & ,nRBpolcomp
     & ,R
     & ,iRlo0,iRlo1
     & ,iRhi0,iRhi1
     & ,B
     & ,iBlo0,iBlo1
     & ,iBhi0,iBhi1
     & ,nBcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      REAL*8 RBtor
      integer nRBpolcomp
      integer iRBpollo0,iRBpollo1
      integer iRBpolhi0,iRBpolhi1
      REAL*8 RBpol(
     & iRBpollo0:iRBpolhi0,
     & iRBpollo1:iRBpolhi1,
     & 0:nRBpolcomp-1)
      integer iRlo0,iRlo1
      integer iRhi0,iRhi1
      REAL*8 R(
     & iRlo0:iRhi0,
     & iRlo1:iRhi1)
      integer nBcomp
      integer iBlo0,iBlo1
      integer iBhi0,iBhi1
      REAL*8 B(
     & iBlo0:iBhi0,
     & iBlo1:iBhi1,
     & 0:nBcomp-1)
      integer i,j
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         B(i,j,0) = RBpol(i,j,0) / R(i,j)
         B(i,j,1) = RBtor / R(i,j)
         B(i,j,2) = RBpol(i,j,1) / R(i,j)
      enddo
      enddo
      return
      end
      subroutine GET_FIELD_MAGNITUDE_AND_UNITVECTOR_2D(
     & igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,Bvec
     & ,iBveclo0,iBveclo1
     & ,iBvechi0,iBvechi1
     & ,nBveccomp
     & ,Bmag
     & ,iBmaglo0,iBmaglo1
     & ,iBmaghi0,iBmaghi1
     & ,bunit
     & ,ibunitlo0,ibunitlo1
     & ,ibunithi0,ibunithi1
     & ,nbunitcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer nBveccomp
      integer iBveclo0,iBveclo1
      integer iBvechi0,iBvechi1
      REAL*8 Bvec(
     & iBveclo0:iBvechi0,
     & iBveclo1:iBvechi1,
     & 0:nBveccomp-1)
      integer iBmaglo0,iBmaglo1
      integer iBmaghi0,iBmaghi1
      REAL*8 Bmag(
     & iBmaglo0:iBmaghi0,
     & iBmaglo1:iBmaghi1)
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1
      integer ibunithi0,ibunithi1
      REAL*8 bunit(
     & ibunitlo0:ibunithi0,
     & ibunitlo1:ibunithi1,
     & 0:nbunitcomp-1)
      integer i,j, l
      double precision sum
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         sum = (0.0d0)
         do l = 0, 2
            sum = sum + Bvec(i,j,l)**2
         enddo
         Bmag(i,j) = sqrt(sum)
         do l = 0, 2
            bunit(i,j,l) = Bvec(i,j,l) / Bmag(i,j)
         enddo
      enddo
      enddo
      return
      end
      subroutine GET_FIELD_DERIVATIVE_DATA_2D(
     & igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,axisymmetric
     & ,RB
     & ,iRBlo0,iRBlo1
     & ,iRBhi0,iRBhi1
     & ,dRBdr
     & ,idRBdrlo0,idRBdrlo1
     & ,idRBdrhi0,idRBdrhi1
     & ,dRBdt
     & ,idRBdtlo0,idRBdtlo1
     & ,idRBdthi0,idRBdthi1
     & ,bunit
     & ,ibunitlo0,ibunitlo1
     & ,ibunithi0,ibunithi1
     & ,nbunitcomp
     & ,dbunitRdr
     & ,idbunitRdrlo0,idbunitRdrlo1
     & ,idbunitRdrhi0,idbunitRdrhi1
     & ,dbunitRdt
     & ,idbunitRdtlo0,idbunitRdtlo1
     & ,idbunitRdthi0,idbunitRdthi1
     & ,dbunitphidr
     & ,idbunitphidrlo0,idbunitphidrlo1
     & ,idbunitphidrhi0,idbunitphidrhi1
     & ,dbunitphidt
     & ,idbunitphidtlo0,idbunitphidtlo1
     & ,idbunitphidthi0,idbunitphidthi1
     & ,dbunitZdr
     & ,idbunitZdrlo0,idbunitZdrlo1
     & ,idbunitZdrhi0,idbunitZdrhi1
     & ,dbunitZdt
     & ,idbunitZdtlo0,idbunitZdtlo1
     & ,idbunitZdthi0,idbunitZdthi1
     & ,R
     & ,iRlo0,iRlo1
     & ,iRhi0,iRhi1
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
     & ,gradb
     & ,igradblo0,igradblo1
     & ,igradbhi0,igradbhi1
     & ,ngradbcomp
     & ,curlb
     & ,icurlblo0,icurlblo1
     & ,icurlbhi0,icurlbhi1
     & ,ncurlbcomp
     & ,bdotcurlb
     & ,ibdotcurlblo0,ibdotcurlblo1
     & ,ibdotcurlbhi0,ibdotcurlbhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer axisymmetric
      integer iRBlo0,iRBlo1
      integer iRBhi0,iRBhi1
      REAL*8 RB(
     & iRBlo0:iRBhi0,
     & iRBlo1:iRBhi1)
      integer idRBdrlo0,idRBdrlo1
      integer idRBdrhi0,idRBdrhi1
      REAL*8 dRBdr(
     & idRBdrlo0:idRBdrhi0,
     & idRBdrlo1:idRBdrhi1)
      integer idRBdtlo0,idRBdtlo1
      integer idRBdthi0,idRBdthi1
      REAL*8 dRBdt(
     & idRBdtlo0:idRBdthi0,
     & idRBdtlo1:idRBdthi1)
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1
      integer ibunithi0,ibunithi1
      REAL*8 bunit(
     & ibunitlo0:ibunithi0,
     & ibunitlo1:ibunithi1,
     & 0:nbunitcomp-1)
      integer idbunitRdrlo0,idbunitRdrlo1
      integer idbunitRdrhi0,idbunitRdrhi1
      REAL*8 dbunitRdr(
     & idbunitRdrlo0:idbunitRdrhi0,
     & idbunitRdrlo1:idbunitRdrhi1)
      integer idbunitRdtlo0,idbunitRdtlo1
      integer idbunitRdthi0,idbunitRdthi1
      REAL*8 dbunitRdt(
     & idbunitRdtlo0:idbunitRdthi0,
     & idbunitRdtlo1:idbunitRdthi1)
      integer idbunitphidrlo0,idbunitphidrlo1
      integer idbunitphidrhi0,idbunitphidrhi1
      REAL*8 dbunitphidr(
     & idbunitphidrlo0:idbunitphidrhi0,
     & idbunitphidrlo1:idbunitphidrhi1)
      integer idbunitphidtlo0,idbunitphidtlo1
      integer idbunitphidthi0,idbunitphidthi1
      REAL*8 dbunitphidt(
     & idbunitphidtlo0:idbunitphidthi0,
     & idbunitphidtlo1:idbunitphidthi1)
      integer idbunitZdrlo0,idbunitZdrlo1
      integer idbunitZdrhi0,idbunitZdrhi1
      REAL*8 dbunitZdr(
     & idbunitZdrlo0:idbunitZdrhi0,
     & idbunitZdrlo1:idbunitZdrhi1)
      integer idbunitZdtlo0,idbunitZdtlo1
      integer idbunitZdthi0,idbunitZdthi1
      REAL*8 dbunitZdt(
     & idbunitZdtlo0:idbunitZdthi0,
     & idbunitZdtlo1:idbunitZdthi1)
      integer iRlo0,iRlo1
      integer iRhi0,iRhi1
      REAL*8 R(
     & iRlo0:iRhi0,
     & iRlo1:iRhi1)
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
      integer ngradbcomp
      integer igradblo0,igradblo1
      integer igradbhi0,igradbhi1
      REAL*8 gradb(
     & igradblo0:igradbhi0,
     & igradblo1:igradbhi1,
     & 0:ngradbcomp-1)
      integer ncurlbcomp
      integer icurlblo0,icurlblo1
      integer icurlbhi0,icurlbhi1
      REAL*8 curlb(
     & icurlblo0:icurlbhi0,
     & icurlblo1:icurlbhi1,
     & 0:ncurlbcomp-1)
      integer ibdotcurlblo0,ibdotcurlblo1
      integer ibdotcurlbhi0,ibdotcurlbhi1
      REAL*8 bdotcurlb(
     & ibdotcurlblo0:ibdotcurlbhi0,
     & ibdotcurlblo1:ibdotcurlbhi1)
      integer i,j, l
      double precision jac, fac
      if (axisymmetric .ne. 0) then
         fac = (1.0d0)
      else
         fac = (0.0d0)
      endif
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         jac = Rr(i,j) * Zt(i,j) - Zr(i,j) * Rt(i,j)
         gradb(i,j,0) = (- RB(i,j)/R(i,j)
     & + (Zt(i,j)*dRBdr(i,j)
     & - Zr(i,j)*dRBdt(i,j)) / jac ) / R(i,j)
         gradb(i,j,1) = (0.0d0)
         gradb(i,j,2) = (-Rt(i,j)*dRBdr(i,j)
     & + Rr(i,j)*dRBdt(i,j)) / (jac * R(i,j))
         curlb(i,j,0) = (Rt(i,j) * dbunitphidr(i,j)
     & - Rr(i,j) * dbunitphidt(i,j)) / jac
         curlb(i,j,1) = (-Rt(i,j) * dbunitRdr(i,j)
     & + Rr(i,j) * dbunitRdt(i,j)
     & - Zt(i,j) * dbunitZdr(i,j)
     & + Zr(i,j) * dbunitZdt(i,j)) / jac
         curlb(i,j,2) = fac * bunit(i,j,1)/R(i,j)
     & + (Zt(i,j) * dbunitphidr(i,j)
     & - Zr(i,j) * dbunitphidt(i,j)) / jac
         bdotcurlb(i,j) = (0.0d0)
         do l = 0, 2
            bdotcurlb(i,j) = bdotcurlb(i,j)
     & + curlb(i,j,l) * bunit(i,j,l)
         end do
      enddo
      enddo
      return
      end
      subroutine DCT_INTERP_2D(
     & coef
     & ,icoeflo0,icoeflo1
     & ,icoefhi0,icoefhi1
     & ,expansion_order
     & ,deriv1
     & ,deriv2
     & ,fac1
     & ,ifac1hi0
     & ,fac2
     & ,ifac2hi0
     & ,sinfac1
     & ,isinfac1hi0
     & ,cosfac1
     & ,icosfac1hi0
     & ,sinfac2
     & ,isinfac2hi0
     & ,cosfac2
     & ,icosfac2hi0
     & ,lambda
     & ,ilambdahi0
     & ,value
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoeflo0,icoeflo1
      integer icoefhi0,icoefhi1
      REAL*8 coef(
     & icoeflo0:icoefhi0,
     & icoeflo1:icoefhi1)
      integer expansion_order
      integer deriv1
      integer deriv2
      integer ifac1hi0
      REAL*8 fac1(
     & 0:ifac1hi0)
      integer ifac2hi0
      REAL*8 fac2(
     & 0:ifac2hi0)
      integer isinfac1hi0
      REAL*8 sinfac1(
     & 0:isinfac1hi0)
      integer icosfac1hi0
      REAL*8 cosfac1(
     & 0:icosfac1hi0)
      integer isinfac2hi0
      REAL*8 sinfac2(
     & 0:isinfac2hi0)
      integer icosfac2hi0
      REAL*8 cosfac2(
     & 0:icosfac2hi0)
      integer ilambdahi0
      REAL*8 lambda(
     & 0:ilambdahi0)
      REAL*8 value
      integer Nminus1, Mminus1, i,j
      Nminus1 = icoefhi0
      Mminus1 = icoefhi1
      value = (0.0d0)
      if (deriv1 .eq. 0 .and. deriv2 .eq. 0) then
         do j = 0, expansion_order - 1
            do i = 0, expansion_order - 1
               value = value + lambda(i) * lambda(j) * coef(i,j) * cosfa
     &c1(i) * cosfac2(j)
            enddo
         enddo
      else if (deriv1 .eq. 1 .and. deriv2 .eq. 0) then
         do j = 0, expansion_order - 1
            do i = 0, expansion_order - 1
               value = value - lambda(i) * lambda(j) * coef(i,j) * sinfa
     &c1(i) * cosfac2(j) * fac1(i)
            enddo
         enddo
      else if (deriv1 .eq. 0 .and. deriv2 .eq. 1) then
         do j = 0, expansion_order - 1
            do i = 0, expansion_order - 1
               value = value - lambda(i) * lambda(j) * coef(i,j) * cosfa
     &c1(i) * sinfac2(j) * fac2(j)
            enddo
         enddo
      else if (deriv1 .eq. 2 .and. deriv2 .eq. 0) then
         do j = 0, expansion_order - 1
            do i = 0, expansion_order - 1
               value = value - lambda(i) * lambda(j) * coef(i,j) * cosfa
     &c1(i) * cosfac2(j) * fac1(i)**2
            enddo
         enddo
      else if (deriv1 .eq. 0 .and. deriv2 .eq. 2) then
         do j = 0, expansion_order - 1
            do i = 0, expansion_order - 1
               value = value - lambda(i) * lambda(j) * coef(i,j) * cosfa
     &c1(i) * cosfac2(j) * fac2(j)**2
            enddo
         enddo
      else if (deriv1 .eq. 1 .and. deriv2 .eq. 1) then
         do j = 0, expansion_order - 1
            do i = 0, expansion_order - 1
               value = value + lambda(i) * lambda(j) * coef(i,j) * sinfa
     &c1(i) * sinfac2(j) * fac1(i) * fac2(j)
            enddo
         enddo
      endif
      value = value * (2.0d0) / dsqrt((Nminus1+(1.0d0))*(Mminus1+(1.0d0)
     &))
      return
      end
