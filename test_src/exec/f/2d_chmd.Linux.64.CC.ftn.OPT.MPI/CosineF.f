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
      subroutine SET_COSINE_2D(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,nxcomp
     & ,dx
     & ,shift
     & ,constant
     & ,amplitude
     & ,mode
     & ,phase
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1,
     & 0:nxcomp-1)
      REAL*8 dx(0:1)
      REAL*8 shift
      REAL*8 constant
      REAL*8 amplitude
      REAL*8 mode(0:1)
      REAL*8 phase(0:1)
      integer i,j
      integer d,n
      REAL*8 arg,product
      do n=0,nphicomp-1
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        x(i,j,0) = (i + (0.500d0))*dx(0) + shift
        product = (1.0d0)
        do d=0,nxcomp-1
          arg = mode(d) * x(i,j,d) + phase(d)
          product = product * cos( arg )
        enddo
        phi(i,j,n) = constant + amplitude * product
      enddo
      enddo
      enddo
      return
      end
