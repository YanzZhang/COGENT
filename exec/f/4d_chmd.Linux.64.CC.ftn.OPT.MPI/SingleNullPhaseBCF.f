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
      subroutine SET_LOGICAL_SHEATH_BC_4D(
     & f
     & ,iflo0,iflo1,iflo2,iflo3
     & ,ifhi0,ifhi1,ifhi2,ifhi3
     & ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2,ibdryboxlo3
     & ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2,ibdryboxhi3
     & ,f_rflct
     & ,if_rflctlo0,if_rflctlo1,if_rflctlo2,if_rflctlo3
     & ,if_rflcthi0,if_rflcthi1,if_rflcthi2,if_rflcthi3
     & ,vel
     & ,ivello0,ivello1,ivello2,ivello3
     & ,ivelhi0,ivelhi1,ivelhi2,ivelhi3
     & ,nvelcomp
     & ,vn
     & ,ivnlo0,ivnlo1,ivnlo2,ivnlo3
     & ,ivnhi0,ivnhi1,ivnhi2,ivnhi3
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,mass
     & ,charge
     & ,iside
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & iflo2:ifhi2,
     & iflo3:ifhi3)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2,ibdryboxlo3
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2,ibdryboxhi3
      integer if_rflctlo0,if_rflctlo1,if_rflctlo2,if_rflctlo3
      integer if_rflcthi0,if_rflcthi1,if_rflcthi2,if_rflcthi3
      REAL*8 f_rflct(
     & if_rflctlo0:if_rflcthi0,
     & if_rflctlo1:if_rflcthi1,
     & if_rflctlo2:if_rflcthi2,
     & if_rflctlo3:if_rflcthi3)
      integer nvelcomp
      integer ivello0,ivello1,ivello2,ivello3
      integer ivelhi0,ivelhi1,ivelhi2,ivelhi3
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & ivello2:ivelhi2,
     & ivello3:ivelhi3,
     & 0:nvelcomp-1)
      integer ivnlo0,ivnlo1,ivnlo2,ivnlo3
      integer ivnhi0,ivnhi1,ivnhi2,ivnhi3
      REAL*8 vn(
     & ivnlo0:ivnhi0,
     & ivnlo1:ivnhi1,
     & ivnlo2:ivnhi2,
     & ivnlo3:ivnhi3)
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3)
      REAL*8 mass
      REAL*8 charge
      integer iside
      integer i,j,k,l,m
      integer isign
      integer vn_jbdry, phi_jbdry
      integer jbdry,jsrc,jsrc_offset
      real pot_energy_min
      if (iside.eq.0) then
        jbdry = ibdryboxhi1
      else
        jbdry = ibdryboxlo1
      endif
      isign = 2*iside-1
      jsrc_offset = 2 * jbdry - isign
      do l = ibdryboxlo3,ibdryboxhi3
      do k = ibdryboxlo2,ibdryboxhi2
      do j = ibdryboxlo1,ibdryboxhi1
      do i = ibdryboxlo0,ibdryboxhi0
          vn_jbdry = jbdry - isign + iside
          if (isign*vn(i,vn_jbdry,k,l).le.(0.0d0)) then
             phi_jbdry = jbdry - isign
             pot_energy_min = -charge * phi(i,phi_jbdry,iphilo2,iphilo3)
             if (mass * (vel(i,jbdry,k,l,0)**2).gt.pot_energy_min) then
               f(i,j,k,l) = (0.0d0)
             else
               jsrc = jsrc_offset - j
               f(i,j,k,l) = f_rflct(i,jsrc,-k-1,l)
             endif
          endif
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine SET_POLOIDAL_DIVERTER_BC_4D(
     & f
     & ,iflo0,iflo1,iflo2,iflo3
     & ,ifhi0,ifhi1,ifhi2,ifhi3
     & ,nfcomp
     & ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2,ibdryboxlo3
     & ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2,ibdryboxhi3
     & ,f_rflct
     & ,if_rflctlo0,if_rflctlo1,if_rflctlo2,if_rflctlo3
     & ,if_rflcthi0,if_rflcthi1,if_rflcthi2,if_rflcthi3
     & ,nf_rflctcomp
     & ,vn
     & ,ivnlo0,ivnlo1,ivnlo2,ivnlo3
     & ,ivnhi0,ivnhi1,ivnhi2,ivnhi3
     & ,nvncomp
     & ,qphi
     & ,iqphilo0,iqphilo1,iqphilo2,iqphilo3
     & ,iqphihi0,iqphihi1,iqphihi2,iqphihi3
     & ,iside
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
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2,ibdryboxlo3
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2,ibdryboxhi3
      integer nf_rflctcomp
      integer if_rflctlo0,if_rflctlo1,if_rflctlo2,if_rflctlo3
      integer if_rflcthi0,if_rflcthi1,if_rflcthi2,if_rflcthi3
      REAL*8 f_rflct(
     & if_rflctlo0:if_rflcthi0,
     & if_rflctlo1:if_rflcthi1,
     & if_rflctlo2:if_rflcthi2,
     & if_rflctlo3:if_rflcthi3,
     & 0:nf_rflctcomp-1)
      integer nvncomp
      integer ivnlo0,ivnlo1,ivnlo2,ivnlo3
      integer ivnhi0,ivnhi1,ivnhi2,ivnhi3
      REAL*8 vn(
     & ivnlo0:ivnhi0,
     & ivnlo1:ivnhi1,
     & ivnlo2:ivnhi2,
     & ivnlo3:ivnhi3,
     & 0:nvncomp-1)
      integer iqphilo0,iqphilo1,iqphilo2,iqphilo3
      integer iqphihi0,iqphihi1,iqphihi2,iqphihi3
      REAL*8 qphi(
     & iqphilo0:iqphihi0,
     & iqphilo1:iqphihi1,
     & iqphilo2:iqphihi2,
     & iqphilo3:iqphihi3)
      integer iside
      integer i,j,k,l,m
      integer isign
      integer n
      integer jbdry,jsrc,jsrc_offset
      REAL*8 pot_energy_min
      if (iside.eq.0) then
        jbdry = ibdryboxhi1
      else
        jbdry = ibdryboxlo1
      endif
      isign = 2*iside-1
      jsrc_offset = 2 * jbdry + isign
      do n=0,nfcomp-1
      do l = ibdryboxlo3,ibdryboxhi3
      do k = ibdryboxlo2,ibdryboxhi2
      do j = ibdryboxlo1,ibdryboxhi1
      do i = ibdryboxlo0,ibdryboxhi0
          pot_energy_min = (-qphi(i,jbdry,iqphilo2,iqphilo3))
          if (isign*vn(i,jbdry,k,l,0).le.(0.0d0)) then
            if ((vn(i,jbdry,k,l,0)**2).gt.pot_energy_min) then
              f(i,j,k,l,n) = (0.0d0)
            else
              jsrc = jsrc_offset - j
              f(i,j,k,l,n) = f_rflct(i,jsrc,-k,l,n)
            endif
          endif
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
