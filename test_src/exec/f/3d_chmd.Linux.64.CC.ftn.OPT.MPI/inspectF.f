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
      subroutine INSPECT_3D(
     & igridboxlo0,igridboxlo1,igridboxlo2
     & ,igridboxhi0,igridboxhi1,igridboxhi2
     & ,data
     & ,idatalo0,idatalo1,idatalo2
     & ,idatahi0,idatahi1,idatahi2
     & ,ndatacomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2
      integer igridboxhi0,igridboxhi1,igridboxhi2
      integer ndatacomp
      integer idatalo0,idatalo1,idatalo2
      integer idatahi0,idatahi1,idatahi2
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1,
     & idatalo2:idatahi2,
     & 0:ndatacomp-1)
      integer i,j,k,n
      do n=0,ndatacomp-1
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
      data(i,j,k,n)=data(i,j,k,n)
      enddo
      enddo
      enddo
      enddo
      return
      end
