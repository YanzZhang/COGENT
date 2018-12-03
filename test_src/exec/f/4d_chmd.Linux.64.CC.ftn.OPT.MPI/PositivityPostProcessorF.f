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
      subroutine REDISTRIBUTENEGATIVES_4D(
     & phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,nphicomp
     & ,iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
     & ,iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
     & ,ineighborBoxlo0,ineighborBoxlo1,ineighborBoxlo2,ineighborBoxlo3
     & ,ineighborBoxhi0,ineighborBoxhi1,ineighborBoxhi2,ineighborBoxhi3
     & ,refVal
     & ,unable
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3,
     & 0:nphicomp-1)
      integer iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorB
     &oxlo3
      integer iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorB
     &oxhi3
      integer ineighborBoxlo0,ineighborBoxlo1,ineighborBoxlo2,ineighborB
     &oxlo3
      integer ineighborBoxhi0,ineighborBoxhi1,ineighborBoxhi2,ineighborB
     &oxhi3
      REAL*8 refVal
      integer unable
      integer n
      integer i,j,k,l
      integer m1,m2,m3,m4
      REAL*8 deltaPhi
      REAL*8 xisum,scale
      REAL*8 xi(-3:3,-3:3,-3:3,-3:3)
      do n=0, (nphicomp-1)
      do l = iinteriorBoxlo3,iinteriorBoxhi3
      do k = iinteriorBoxlo2,iinteriorBoxhi2
      do j = iinteriorBoxlo1,iinteriorBoxhi1
      do i = iinteriorBoxlo0,iinteriorBoxhi0
          if (phi(i,j,k,l,n).lt.(0.0d0)) then
            deltaPhi = -phi(i,j,k,l,n)
            if ((refVal+deltaPhi).ne.refVal) then
              xisum = (0.0d0)
      do m4 = ineighborBoxlo3,ineighborBoxhi3
      do m3 = ineighborBoxlo2,ineighborBoxhi2
      do m2 = ineighborBoxlo1,ineighborBoxhi1
      do m1 = ineighborBoxlo0,ineighborBoxhi0
                 xi(m1,m2,m3,m4)
     & = max((0.0d0),phi(i+m1,j+m2,k+m3,l+m4,n))
                 xisum = xisum + xi(m1,m2,m3,m4)
      enddo
      enddo
      enddo
      enddo
              if (xisum.ge.deltaPhi) then
                scale = deltaPhi / xisum
      do m4 = ineighborBoxlo3,ineighborBoxhi3
      do m3 = ineighborBoxlo2,ineighborBoxhi2
      do m2 = ineighborBoxlo1,ineighborBoxhi1
      do m1 = ineighborBoxlo0,ineighborBoxhi0
                   xi(m1,m2,m3,m4)
     & = xi(m1,m2,m3,m4) * scale
      enddo
      enddo
      enddo
      enddo
      do m4 = ineighborBoxlo3,ineighborBoxhi3
      do m3 = ineighborBoxlo2,ineighborBoxhi2
      do m2 = ineighborBoxlo1,ineighborBoxhi1
      do m1 = ineighborBoxlo0,ineighborBoxhi0
                   phi(i+m1,j+m2,k+m3,l+m4,n)
     & = phi(i+m1,j+m2,k+m3,l+m4,n)
     & - xi(m1,m2,m3,m4)
      enddo
      enddo
      enddo
      enddo
                phi(i,j,k,l,n) = (0.0d0)
              else
                unable = unable + 1
              endif
            else
              phi(i,j,k,l,n) = (0.0d0)
            endif
          endif
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FINDANYNEGATIVES_4D(
     & count
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,nphicomp
     & ,iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
     & ,iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer count
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3,
     & 0:nphicomp-1)
      integer iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorB
     &oxlo3
      integer iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorB
     &oxhi3
      integer n
      integer i,j,k,l
      REAL*8 val
      do n=0, (nphicomp-1)
      do l = iinteriorBoxlo3,iinteriorBoxhi3
      do k = iinteriorBoxlo2,iinteriorBoxhi2
      do j = iinteriorBoxlo1,iinteriorBoxhi1
      do i = iinteriorBoxlo0,iinteriorBoxhi0
          if (phi(i,j,k,l,n).lt.(0.0d0)) then
            count = count + 1
            return
          endif
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FINDALLNEGATIVES_4D(
     & count
     & ,phi
     & ,iphilo0,iphilo1,iphilo2,iphilo3
     & ,iphihi0,iphihi1,iphihi2,iphihi3
     & ,nphicomp
     & ,iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
     & ,iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer count
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & iphilo3:iphihi3,
     & 0:nphicomp-1)
      integer iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorB
     &oxlo3
      integer iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorB
     &oxhi3
      integer n
      integer i,j,k,l
      do n=0, (nphicomp-1)
      do l = iinteriorBoxlo3,iinteriorBoxhi3
      do k = iinteriorBoxlo2,iinteriorBoxhi2
      do j = iinteriorBoxlo1,iinteriorBoxhi1
      do i = iinteriorBoxlo0,iinteriorBoxhi0
          if (phi(i,j,k,l,n).lt.(0.0d0)) then
            count = count + 1
          endif
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
