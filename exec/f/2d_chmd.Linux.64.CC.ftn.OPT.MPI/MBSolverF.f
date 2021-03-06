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
      subroutine ACCUM_FLUX_STENCIL4_2D(
     & dir
     & ,deriv_dir
     & ,side
     & ,h
     & ,coef
     & ,icoeflo0,icoeflo1
     & ,icoefhi0,icoefhi1
     & ,global
     & ,sum
     & ,isumlo0,isumlo1
     & ,isumhi0,isumhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer dir
      integer deriv_dir
      integer side
      REAL*8 h(0:1)
      integer icoeflo0,icoeflo1
      integer icoefhi0,icoefhi1
      REAL*8 coef(
     & icoeflo0:icoefhi0,
     & icoeflo1:icoefhi1)
      integer global(0:1)
      integer isumlo0,isumlo1
      integer isumhi0,isumhi1
      REAL*8 sum(
     & isumlo0:isumhi0,
     & isumlo1:isumhi1)
      integer i0,i1, ii0,ii1, it0,it1, ic0,ic1, lo0,lo1,
     & hi0,hi1, nlo0,nlo1, nhi0,nhi1, tii0,tii1,
     & tdir, m, n, n_start, n_stop, t_start, t_stop
      REAL*8 hfac, n_stencil(0:3), t_stencil_1(0:2 -1,0:4), t_stencil_2(
     &0:2 -1,0:4),
     & d, trans_grad_d, t_stencil_1_prod, t_stencil_2_prod
      ic0 = (isumlo0+isumhi0)/2
      ic1 = (isumlo1+isumhi1)/2
      hfac = dble(1 - 2*side) / h(deriv_dir)
      do tdir = 0, 2 -1
         if (tdir .ne. dir) then
            hfac = hfac * h(tdir)
         endif
      enddo
      ii0=CHF_ID(0,dir)
      ii1=CHF_ID(1,dir)
      d = coef(global(0)+side*ii0,global(1)+side*ii1)
      if (dir .eq. deriv_dir) then
         n_stencil(0) = (1.0d0) / 24.d0
         n_stencil(1) = -27.d0 / 24.d0
         n_stencil(2) = 27.d0 / 24.d0
         n_stencil(3) = -(1.0d0) / 24.d0
      else
         n_stencil(0) = -(1.0d0) / 16.d0
         n_stencil(1) = (9.0d0) / 16.d0
         n_stencil(2) = (9.0d0) / 16.d0
         n_stencil(3) = -(1.0d0) / 16.d0
      endif
      n_start = side - 2
      n_stop = n_start + 3
      nlo0 = ic0 + ii0*n_start
      nhi0 = ic0 + ii0*n_stop
      nlo1 = ic1 + ii1*n_start
      nhi1 = ic1 + ii1*n_stop
      do m = 0, 4
         do n = 0, 2 -1
            t_stencil_1(n,m) = (0.0d0)
            t_stencil_2(n,m) = (0.0d0)
         enddo
      enddo
      if (dir .eq. deriv_dir) then
         lo0 = nlo0
         hi0 = nhi0
         lo1 = nlo1
         hi1 = nhi1
            do i1 = lo1, hi1
         do i0 = lo0, hi0
                  n = ii0*(i0-lo0) + ii1*(i1-lo1)
                  sum(i0,i1) = sum(i0,i1) + hfac * d * n_stencil(n)
            enddo
               enddo
         t_start = -1
         t_stop = 1
         do tdir = 0, 2 -1
            if (tdir .ne. dir) then
               do n = 0, 2 -1
                  t_stencil_2(n,0) = (0.0d0)
               enddo
               t_stencil_2(tdir,0) = (1.0d0)
               t_stencil_2(tdir,1) = -(2.0d0)
               t_stencil_2(tdir,2) = (1.0d0)
      tii0=CHF_ID(0,tdir)
      tii1=CHF_ID(1,tdir)
               lo0 = nlo0 + tii0*t_start
               hi0 = nhi0 + tii0*t_stop
               lo1 = nlo1 + tii1*t_start
               hi1 = nhi1 + tii1*t_stop
                  do i1 = lo1, hi1
               do i0 = lo0, hi0
                        n = ii0*(i0-lo0) + ii1*(i1-lo1)
                        t_stencil_2_prod = 0 + ii0*t_stencil_2(1,i1-lo1)
     & + ii1*t_stencil_2(0,i0-lo0)
                        sum(i0,i1) = sum(i0,i1)
     & + hfac * d * n_stencil(n) * t_stencil_2_prod / 24.d0
                  enddo
                     enddo
            endif
         enddo
      else
         do n = 0, 2 -1
            t_stencil_1(n,0) = (0.0d0)
            t_stencil_2(n,0) = (0.0d0)
         enddo
         t_stencil_1(deriv_dir,0) = (1.0d0) / (12.0d0)
         t_stencil_1(deriv_dir,1) = -(8.0d0) / (12.0d0)
         t_stencil_1(deriv_dir,2) = (0.0d0)
         t_stencil_1(deriv_dir,3) = (8.0d0) / (12.0d0)
         t_stencil_1(deriv_dir,4) = -(1.0d0) / (12.0d0)
         t_stencil_2(deriv_dir,0) = -(0.500d0)
         t_stencil_2(deriv_dir,1) = (1.0d0)
         t_stencil_2(deriv_dir,2) = (0.0d0)
         t_stencil_2(deriv_dir,3) = -(1.0d0)
         t_stencil_2(deriv_dir,4) = (0.500d0)
         t_start = -2
         t_stop = 2
      tii0=CHF_ID(0,deriv_dir)
      tii1=CHF_ID(1,deriv_dir)
         lo0 = nlo0 + tii0*t_start
         hi0 = nhi0 + tii0*t_stop
         lo1 = nlo1 + tii1*t_start
         hi1 = nhi1 + tii1*t_stop
            do i1 = lo1, hi1
         do i0 = lo0, hi0
                  n = ii0*(i0-lo0) + ii1*(i1-lo1)
                  t_stencil_1_prod = 0 + ii0*t_stencil_1(1,i1-lo1) + ii1
     &*t_stencil_1(0,i0-lo0)
                  t_stencil_2_prod = 0 + ii0*t_stencil_2(1,i1-lo1) + ii1
     &*t_stencil_2(0,i0-lo0)
                  sum(i0,i1) = sum(i0,i1)
     & + hfac * d * n_stencil(n) * (t_stencil_1_prod + t_stencil_2_prod 
     &/ 24.d0)
            enddo
               enddo
      endif
      if (dir .eq. deriv_dir) then
         n_stencil(0) = -(1.0d0)
         n_stencil(1) = (1.0d0)
      else
         n_stencil(0) = (0.500d0)
         n_stencil(1) = (0.500d0)
      endif
      n_start = side - 1
      n_stop = n_start + 1
      nlo0 = ic0 + ii0*n_start
      nhi0 = ic0 + ii0*n_stop
      nlo1 = ic1 + ii1*n_start
      nhi1 = ic1 + ii1*n_stop
      t_start = -1
      t_stop = 1
      do tdir = 0, 2 -1
         if (tdir .ne. dir) then
            do n = 0, 2 -1
               t_stencil_2(n,0) = (0.0d0)
            enddo
            t_stencil_2(tdir,0) = -(0.500d0)
            t_stencil_2(tdir,1) = (0.0d0)
            t_stencil_2(tdir,2) = (0.500d0)
      tii0=CHF_ID(0,tdir)
      tii1=CHF_ID(1,tdir)
            lo0 = nlo0 + tii0*t_start
            hi0 = nhi0 + tii0*t_stop
            lo1 = nlo1 + tii1*t_start
            hi1 = nhi1 + tii1*t_stop
      it0=CHF_ID(0,tdir)
      it1=CHF_ID(1,tdir)
            trans_grad_d = (0.500d0) * (
     & coef(global(0)+ii0*side+it0,global(1)+ii1*side+it1)
     & - coef(global(0)+ii0*side-it0,global(1)+ii1*side-it1)
     & )
               do i1 = lo1, hi1
            do i0 = lo0, hi0
                     n = ii0*(i0-lo0) + ii1*(i1-lo1)
                     t_stencil_2_prod = 0 + ii0*t_stencil_2(1,i1-lo1) + 
     &ii1*t_stencil_2(0,i0-lo0)
                     sum(i0,i1) = sum(i0,i1)
     & + hfac * trans_grad_d * n_stencil(n) * t_stencil_2_prod / (12.0d0
     &)
               enddo
                  enddo
         endif
      enddo
      return
      end
      subroutine ACCUM_FLUX_STENCIL2_2D(
     & dir
     & ,deriv_dir
     & ,side
     & ,h
     & ,coef
     & ,icoeflo0,icoeflo1
     & ,icoefhi0,icoefhi1
     & ,global
     & ,sum
     & ,isumlo0,isumlo1
     & ,isumhi0,isumhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer dir
      integer deriv_dir
      integer side
      REAL*8 h(0:1)
      integer icoeflo0,icoeflo1
      integer icoefhi0,icoefhi1
      REAL*8 coef(
     & icoeflo0:icoefhi0,
     & icoeflo1:icoefhi1)
      integer global(0:1)
      integer isumlo0,isumlo1
      integer isumhi0,isumhi1
      REAL*8 sum(
     & isumlo0:isumhi0,
     & isumlo1:isumhi1)
      integer i0,i1, ii0,ii1, ic0,ic1, tii0,tii1,
     & lo0,lo1, hi0,hi1, nlo0,nlo1, nhi0,nhi1,
     & tdir, m, n, n_start, n_stop, t_start, t_stop
      REAL*8 hfac, d, n_stencil(0:1), t_stencil(0:2 -1,0:2), t_stencil_p
     &rod
      ic0 = (isumlo0+isumhi0)/2
      ic1 = (isumlo1+isumhi1)/2
      hfac = (1 - 2*side)/ h(deriv_dir)
      do tdir = 0, 2 -1
         if (tdir .ne. dir) then
            hfac = hfac * h(tdir)
         endif
      enddo
      ii0=CHF_ID(0,dir)
      ii1=CHF_ID(1,dir)
      d = coef(global(0)+side*ii0,global(1)+side*ii1)
      if (dir .eq. deriv_dir) then
         n_stencil(0) = -(1.0d0)
         n_stencil(1) = (1.0d0)
      else
         n_stencil(0) = (0.500d0)
         n_stencil(1) = (0.500d0)
      endif
      n_start = side - 1
      n_stop = n_start + 1
      nlo0 = ic0 + ii0*n_start
      nhi0 = ic0 + ii0*n_stop
      nlo1 = ic1 + ii1*n_start
      nhi1 = ic1 + ii1*n_stop
      do m = 0, 2
         do n = 0, 2 -1
            t_stencil(n,m) = (0.0d0)
         enddo
      enddo
      if (dir .eq. deriv_dir) then
         lo0 = nlo0
         hi0 = nhi0
         lo1 = nlo1
         hi1 = nhi1
            do i1 = lo1, hi1
         do i0 = lo0, hi0
                  n = ii0*(i0-lo0) + ii1*(i1-lo1)
                  sum(i0,i1) = sum(i0,i1) + hfac * d * n_stencil(n)
            enddo
               enddo
      else
         do n = 0, 2 -1
            t_stencil(n,0) = (0.0d0)
         enddo
         t_stencil(deriv_dir,0) = -(0.500d0)
         t_stencil(deriv_dir,1) = (0.0d0)
         t_stencil(deriv_dir,2) = (0.500d0)
         t_start = -1
         t_stop = 1
      tii0=CHF_ID(0,deriv_dir)
      tii1=CHF_ID(1,deriv_dir)
         lo0 = nlo0 + tii0*t_start
         hi0 = nhi0 + tii0*t_stop
         lo1 = nlo1 + tii1*t_start
         hi1 = nhi1 + tii1*t_stop
            do i1 = lo1, hi1
         do i0 = lo0, hi0
                  n = ii0*(i0-lo0) + ii1*(i1-lo1)
                  t_stencil_prod = 0 + ii0*t_stencil(1,i1-lo1) + ii1*t_s
     &tencil(0,i0-lo0)
                  sum(i0,i1) = sum(i0,i1) + hfac * d * n_stencil(n) * t_
     &stencil_prod
            enddo
               enddo
      endif
      return
      end
