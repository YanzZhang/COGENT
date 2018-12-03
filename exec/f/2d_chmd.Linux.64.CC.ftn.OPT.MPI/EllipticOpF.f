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
      subroutine CELL_CENTERED_FIELD_COMPONENT_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dir
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,h
     & ,order
     & ,Efield
     & ,iEfieldlo0,iEfieldlo1
     & ,iEfieldhi0,iEfieldhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      REAL*8 h(0:1)
      integer order
      integer iEfieldlo0,iEfieldlo1
      integer iEfieldhi0,iEfieldhi1
      REAL*8 Efield(
     & iEfieldlo0:iEfieldhi0,
     & iEfieldlo1:iEfieldhi1)
      integer i,j, ii,jj
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      if (order .eq. 4) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            Efield(i,j) = -(
     & phi(i-2*ii,j-2*jj)
     & - (8.0d0) * phi(i- ii,j- jj)
     & + (8.0d0) * phi(i+ ii,j+ jj)
     & - phi(i+2*ii,j+2*jj)
     & ) / ((12.0d0) * h(dir))
      enddo
      enddo
      else if (order .eq. 2) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            Efield(i,j) = -(
     & - phi(i- ii,j- jj)
     & + phi(i+ ii,j+ jj)
     & ) / ((2.0d0) * h(dir))
      enddo
      enddo
      endif
      return
      end
      subroutine FACE_CENTERED_FIELD_COMPONENT_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dir
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,h
     & ,order
     & ,Efield
     & ,iEfieldlo0,iEfieldlo1
     & ,iEfieldhi0,iEfieldhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      REAL*8 h(0:1)
      integer order
      integer iEfieldlo0,iEfieldlo1
      integer iEfieldhi0,iEfieldhi1
      REAL*8 Efield(
     & iEfieldlo0:iEfieldhi0,
     & iEfieldlo1:iEfieldhi1)
      integer i,j, ii,jj
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      if (order .eq. 4) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            Efield(i,j) = -(
     & 27.d0 * (phi(i ,j )
     & - phi(i- ii,j- jj))
     & - phi(i+ ii,j+ jj)
     & + phi(i-2*ii,j-2*jj)
     & ) / (24.d0 * h(dir))
      enddo
      enddo
      else if (order .eq. 2) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            Efield(i,j) = -(
     & - phi(i- ii,j- jj)
     & + phi(i ,j )
     & ) / h(dir)
      enddo
      enddo
      endif
      return
      end
      subroutine FACE_INTERPOLATE_2D(
     & dir
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,order
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,data
     & ,idatalo0,idatalo1
     & ,idatahi0,idatahi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer dir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer order
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1)
      integer i,j, ii,jj
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      if (order .eq. 4) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            data(i,j) = (
     & 9.d0 * (phi(i ,j )
     & + phi(i- ii,j- jj))
     & - (phi(i+ ii,j+ jj)
     & + phi(i-2*ii,j-2*jj))
     & ) / 16.d0
      enddo
      enddo
      else if (order .eq. 2) then
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         data(i,j) = (
     & phi(i- ii,j- jj)
     & + phi(i ,j )
     & ) / 2.d0
      enddo
      enddo
      endif
      return
      end
      subroutine EXTRAP_FOR_CC_OPS_2D(
     & dir
     & ,side
     & ,order
     & ,ifaceboxlo0,ifaceboxlo1
     & ,ifaceboxhi0,ifaceboxhi1
     & ,iinteriorboxlo0,iinteriorboxlo1
     & ,iinteriorboxhi0,iinteriorboxhi1
     & ,array
     & ,iarraylo0,iarraylo1
     & ,iarrayhi0,iarrayhi1
     & ,narraycomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer dir
      integer side
      integer order
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer iinteriorboxlo0,iinteriorboxlo1
      integer iinteriorboxhi0,iinteriorboxhi1
      integer narraycomp
      integer iarraylo0,iarraylo1
      integer iarrayhi0,iarrayhi1
      REAL*8 array(
     & iarraylo0:iarrayhi0,
     & iarraylo1:iarrayhi1,
     & 0:narraycomp-1)
      integer i,id,ni,j,jd,nj, m, n, comp, ncomp
      double precision sum, coef2(0:2,2), coef4(0:4,3)
      data coef2 / 3.d0, -3.d0, 1.d0,
     & 6.d0, -8.d0, 3.d0 /
      data coef4 / 5.d0, -10.d0, 10.d0, -5.d0, 1.d0,
     & 15.d0, -40.d0, 45.d0, -24.d0, 5.d0,
     & 35.d0, -105.d0, 126.d0, -70.d0, 15.d0 /
      ncomp = narraycomp
      id = CHF_ID(0,dir)*side
      jd = CHF_ID(1,dir)*side
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
        if (side .eq. -1) then
           ni = id*(i-iinteriorboxlo0)
           n = ni
           nj = jd*(j-iinteriorboxlo1)
           n = n + nj
         else if (side .eq. 1) then
            ni = id*(i-iinteriorboxhi0)
            n = ni
            nj = jd*(j-iinteriorboxhi1)
            n = n + nj
          endif
          do comp = 0, ncomp-1
             sum = (0.0d0)
             if (order .eq. 4) then
                do m = 0, 4
                   sum = sum + coef4(m,n)*array(i-id*(ni+m),j-jd*(nj+m),
     &comp)
                enddo
                array(i,j,comp) = sum
             else if (order .eq. 2) then
                do m = 0, 2
                   sum = sum + coef2(m,n)*array(i-id*(ni+m),j-jd*(nj+m),
     &comp)
                enddo
                array(i,j,comp) = sum
             endif
          enddo
      enddo
      enddo
      return
      end
      subroutine EXTRAP_FOR_FC_OPS_2D(
     & dir
     & ,side
     & ,order
     & ,ifaceboxlo0,ifaceboxlo1
     & ,ifaceboxhi0,ifaceboxhi1
     & ,iinteriorboxlo0,iinteriorboxlo1
     & ,iinteriorboxhi0,iinteriorboxhi1
     & ,array
     & ,iarraylo0,iarraylo1
     & ,iarrayhi0,iarrayhi1
     & ,narraycomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer dir
      integer side
      integer order
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer iinteriorboxlo0,iinteriorboxlo1
      integer iinteriorboxhi0,iinteriorboxhi1
      integer narraycomp
      integer iarraylo0,iarraylo1
      integer iarrayhi0,iarrayhi1
      REAL*8 array(
     & iarraylo0:iarrayhi0,
     & iarraylo1:iarrayhi1,
     & 0:narraycomp-1)
      integer i,id,ni,j,jd,nj, m, n, comp, ncomp
      double precision sum, coef2(0:2,2), coef4(0:4,2)
      data coef2 / 3.d0, -3.d0, 1.d0,
     & 6.d0, -8.d0, 3.d0 /
      data coef4 / 4.625d0, -8.5d0, 7.75d0, -3.5d0, 0.625d0,
     & 11.25d0, -25.d0, 22.5d0, -9.d0, 1.25d0 /
      ncomp = narraycomp
      id = CHF_ID(0,dir)*side
      jd = CHF_ID(1,dir)*side
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
        if (side .eq. -1) then
           ni = id*(i-iinteriorboxlo0)
           n = ni
           nj = jd*(j-iinteriorboxlo1)
           n = n + nj
         else if (side .eq. 1) then
            ni = id*(i-iinteriorboxhi0)
            n = ni
            nj = jd*(j-iinteriorboxhi1)
            n = n + nj
          endif
          do comp = 0, ncomp-1
             sum = (0.0d0)
             if (order .eq. 4) then
                do m = 0, 4
                   sum = sum + coef4(m,n)*array(i-id*(ni+m),j-jd*(nj+m),
     &comp)
                enddo
                array(i,j,comp) = sum
             else if (order .eq. 2) then
                do m = 0, 2
                   sum = sum + coef2(m,n)*array(i-id*(ni+m),j-jd*(nj+m),
     &comp)
                enddo
                array(i,j,comp) = sum
             endif
          enddo
      enddo
      enddo
      return
      end
