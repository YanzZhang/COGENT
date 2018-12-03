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
      subroutine SET_MAXWELL4D_4D(
     & f
     & ,iflo0,iflo1,iflo2,iflo3
     & ,ifhi0,ifhi1,ifhi2,ifhi3
     & ,nfcomp
     & ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     & ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     & ,coords
     & ,icoordslo0,icoordslo1,icoordslo2,icoordslo3
     & ,icoordshi0,icoordshi1,icoordshi2,icoordshi3
     & ,ncoordscomp
     & ,dens
     & ,idenslo0,idenslo1,idenslo2,idenslo3
     & ,idenshi0,idenshi1,idenshi2,idenshi3
     & ,temp
     & ,itemplo0,itemplo1,itemplo2,itemplo3
     & ,itemphi0,itemphi1,itemphi2,itemphi3
     & ,vshift
     & ,ivshiftlo0,ivshiftlo1,ivshiftlo2,ivshiftlo3
     & ,ivshifthi0,ivshifthi1,ivshifthi2,ivshifthi3
     & ,b
     & ,iblo0,iblo1,iblo2,iblo3
     & ,ibhi0,ibhi1,ibhi2,ibhi3
     & ,mass
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfcomp
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & iflo2:ifhi2,
     & iflo3:ifhi3,
     & 0:nfcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer ncoordscomp
      integer icoordslo0,icoordslo1,icoordslo2,icoordslo3
      integer icoordshi0,icoordshi1,icoordshi2,icoordshi3
      REAL*8 coords(
     & icoordslo0:icoordshi0,
     & icoordslo1:icoordshi1,
     & icoordslo2:icoordshi2,
     & icoordslo3:icoordshi3,
     & 0:ncoordscomp-1)
      integer idenslo0,idenslo1,idenslo2,idenslo3
      integer idenshi0,idenshi1,idenshi2,idenshi3
      REAL*8 dens(
     & idenslo0:idenshi0,
     & idenslo1:idenshi1,
     & idenslo2:idenshi2,
     & idenslo3:idenshi3)
      integer itemplo0,itemplo1,itemplo2,itemplo3
      integer itemphi0,itemphi1,itemphi2,itemphi3
      REAL*8 temp(
     & itemplo0:itemphi0,
     & itemplo1:itemphi1,
     & itemplo2:itemphi2,
     & itemplo3:itemphi3)
      integer ivshiftlo0,ivshiftlo1,ivshiftlo2,ivshiftlo3
      integer ivshifthi0,ivshifthi1,ivshifthi2,ivshifthi3
      REAL*8 vshift(
     & ivshiftlo0:ivshifthi0,
     & ivshiftlo1:ivshifthi1,
     & ivshiftlo2:ivshifthi2,
     & ivshiftlo3:ivshifthi3)
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL*8 b(
     & iblo0:ibhi0,
     & iblo1:ibhi1,
     & iblo2:ibhi2,
     & iblo3:ibhi3)
      REAL*8 mass
      integer i,j,k,l,n
      integer ivp,imu
      REAL*8 denloc, temploc, vshiftloc, bloc
      REAL*8 vpar, mu
      REAL*8 eparnorm, munorm
      REAL*8 factor, val
      REAL*8 minf, maxf
      minf=1.0d30
      maxf=(0.0d0)
      factor = dsqrt((3.14159265358979323846264338327950288D0)*((2.0d0)/
     &mass)**3)
      ivp = ncoordscomp-2
      imu = ncoordscomp-1
      do n=0,nfcomp-1
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         vshiftloc = vshift(i,j,itemplo2,itemplo3)
         temploc = temp(i,j,itemplo2,itemplo3)
         denloc = dens(i,j,idenslo2,idenslo3)
         bloc = b(i,j,iblo2,iblo3)
         vpar = coords(i,j,k,l,ivp)
         mu = coords(i,j,k,l,imu)
         eparnorm = (0.500d0) * mass * (vpar-vshiftloc)**2 / temploc
         munorm = (0.500d0) * bloc * mu / temploc
         val = dexp( -( eparnorm + munorm ) )
         val = val * denloc / ( factor * dsqrt( temploc) * temploc )
         minf = min(minf,val)
         maxf = max(maxf,val)
         f(i,j,k,l,n) = val
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
