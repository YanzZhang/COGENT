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
      subroutine SET_POLOIDAL_BC_4D(
     & iboundaryboxlo0,iboundaryboxlo1,iboundaryboxlo2,iboundaryboxlo3
     & ,iboundaryboxhi0,iboundaryboxhi1,iboundaryboxhi2,iboundaryboxhi3
     & ,flipdist
     & ,iflipdistlo0,iflipdistlo1,iflipdistlo2,iflipdistlo3
     & ,iflipdisthi0,iflipdisthi1,iflipdisthi2,iflipdisthi3
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,vp
     & ,ivphi0
     & ,vp_lo
     & ,psi
     & ,ipsihi0
     & ,psilo
     & ,psi_lim
     & ,f
     & ,iflo0,iflo1,iflo2,iflo3
     & ,ifhi0,ifhi1,ifhi2,ifhi3
     & ,charge
     & ,ilohi
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboundaryboxlo0,iboundaryboxlo1,iboundaryboxlo2,iboundaryb
     &oxlo3
      integer iboundaryboxhi0,iboundaryboxhi1,iboundaryboxhi2,iboundaryb
     &oxhi3
      integer iflipdistlo0,iflipdistlo1,iflipdistlo2,iflipdistlo3
      integer iflipdisthi0,iflipdisthi1,iflipdisthi2,iflipdisthi3
      REAL*8 flipdist(
     & iflipdistlo0:iflipdisthi0,
     & iflipdistlo1:iflipdisthi1,
     & iflipdistlo2:iflipdisthi2,
     & iflipdistlo3:iflipdisthi3)
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3)
      integer ivphi0
      REAL*8 vp(
     & 0:ivphi0)
      integer vp_lo
      integer ipsihi0
      REAL*8 psi(
     & 0:ipsihi0)
      integer psilo
      REAL*8 psi_lim
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & iflo2:ifhi2,
     & iflo3:ifhi3)
      REAL*8 charge
      integer ilohi
      double precision vploc,vpsq,df,phibndy
      integer i,j,k,l
      integer jbndy,jbndy1,vp_hi,nvp,twojbndy1,jsrc,jbndy2,jbndygtr,jbnd
     &yless,iparsign
      if (ilohi == 1) then
           jbndy = iboundaryboxhi1
           jbndy1 = jbndy+1
           jbndy2 = jbndy+2
           jbndygtr = jbndy2
           jbndyless = jbndy1
           iparsign = 1
      else
           jbndy = iboundaryboxlo1
           jbndy1 = jbndy-1
           jbndy2 = jbndy-2
           jbndygtr = jbndy1
           jbndyless = jbndy2
           iparsign = -1
      endif
      twojbndy1 = jbndy+jbndy1
      do l = iboundaryboxlo3,iboundaryboxhi3
      do k = iboundaryboxlo2,iboundaryboxhi2
      do j = iboundaryboxlo1,iboundaryboxhi1
      do i = iboundaryboxlo0,iboundaryboxhi0
         phibndy = phi(i,jbndy,iphilo2,iphilo3)
         if (psi(i-psilo) .ge. psi_lim) then
            df = f(i,jbndyless,k,l)-f(i,jbndygtr,k,l)
            vploc = iparsign*vp(k-vp_lo)
            vpsq=vploc**2
            if (1 .lt. 0) then
                f(i,j,k,l) = f(i,jbndy1,k,l)+(jbndy1-j)*df
            else
                jsrc = twojbndy1-j
                vpsq = vploc**2
                if (vpsq .gt. (-charge*phibndy)) then
                     f(i,j,k,l) = 0.
                else
                    f(i,j,k,l) = flipdist(i,jsrc,-k,l)
                endif
            endif
         endif
      enddo
      enddo
      enddo
      enddo
      return
      end
