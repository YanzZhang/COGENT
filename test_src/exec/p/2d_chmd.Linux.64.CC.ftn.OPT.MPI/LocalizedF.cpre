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
      subroutine SET_LOCALIZED_2D(
     & f
     & ,iflo0,iflo1
     & ,ifhi0,ifhi1
     & ,nfcomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,coord
     & ,icoordlo0,icoordlo1
     & ,icoordhi0,icoordhi1
     & ,ncoordcomp
     & ,amp
     & ,center
     & ,width
     & ,floor
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfcomp
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & 0:nfcomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ncoordcomp
      integer icoordlo0,icoordlo1
      integer icoordhi0,icoordhi1
      REAL*8 coord(
     & icoordlo0:icoordhi0,
     & icoordlo1:icoordhi1,
     & 0:ncoordcomp-1)
      REAL*8 amp
      REAL*8 center(0:1)
      REAL*8 width(0:1)
      REAL*8 floor
      integer i,j
      integer d,n
      REAL*8 beta, x
      do n=0,nfcomp-1
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        beta = (0.0d0)
        do d=0,ncoordcomp-1
          x = coord(i,j,d) - center(d)
          beta = beta + ( x / width(d) )**2
        enddo
        f(i,j,n) = amp * dexp(-beta) + floor
      enddo
      enddo
      enddo
      return
      end
