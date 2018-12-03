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
      subroutine GET_SLAB_FIELD_DATA_2D(
     & igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,xix
     & ,ixixlo0,ixixlo1
     & ,ixixhi0,ixixhi1
     & ,ByInner
     & ,ByOuter
     & ,BzInner
     & ,BzOuter
     & ,xmax
     & ,ximax
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
      integer ixixlo0,ixixlo1
      integer ixixhi0,ixixhi1
      REAL*8 xix(
     & ixixlo0:ixixhi0,
     & ixixlo1:ixixhi1)
      REAL*8 ByInner
      REAL*8 ByOuter
      REAL*8 BzInner
      REAL*8 BzOuter
      REAL*8 xmax
      REAL*8 ximax
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
      integer i0,i1, l
      double precision
     & xi, Bz, By, BzBy, dBydx, dBzdx, Bmag, bunitDotcurlb,
     & bunit(0:2), Bvec(0:2), gradB(0:2), curlb(0:2)
      dBydx = (ByOuter - ByInner)/xmax
      dBzdx = (BzOuter - BzInner)/xmax
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
         xi = xix(i0,i1)
         By = ByInner + dBydx * xi * (xmax/ximax)
         Bz = BzInner + dBzdx * xi * (xmax/ximax)
         Bvec(0) = 0.
         Bvec(1) = By
         Bvec(2) = Bz
         Bmag = dsqrt(By*By + Bz*Bz)
         bunit(0) = 0.
         bunit(1) = By/Bmag
         bunit(2) = Bz/Bmag
         gradB(0) = (Bz*dBzdx + By*dBydx)/Bmag
         gradB(1) = (0.0d0)
         gradB(2) = (0.0d0)
         curlb(0) = (0.0d0)
         curlb(1) = -(By**2 * dBzdx - By*Bz*dBydx)/(Bmag**3)
         curlb(2) = (Bz**2 * dBydx - Bz*By*dBzdx)/(Bmag**3)
         if (curlb(1)**2 + curlb(2)**2 .gt. 1.0e-10) then
           bunitDotcurlb = bunit(1) * curlb(1) + bunit(2) * curlb(2)
           bunitDotcurlb = bunitDotcurlb / dsqrt(curlb(1)**2 + curlb(2)*
     &*2)
           bunitDotcurlb = (0.0d0)
         endif
         Bmag_pt(i0,i1) = Bmag
         bdotcurlb_pt(i0,i1) = bunitDotcurlb
         do l = 0, 2
            bunit_pt(i0,i1,l) = bunit(l)
            b_pt(i0,i1,l) = Bvec(l)
            gradb_pt(i0,i1,l) = gradB(l)
            curlb_pt(i0,i1,l) = curlb(l)
         end do
      enddo
      enddo
      return
      end
