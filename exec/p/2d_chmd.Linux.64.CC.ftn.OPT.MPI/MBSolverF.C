#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine ACCUM_FLUX_STENCIL4_2D(
     &           dir
     &           ,deriv_dir
     &           ,side
     &           ,h
     &           ,coef
     &           ,icoeflo0,icoeflo1
     &           ,icoefhi0,icoefhi1
     &           ,global
     &           ,sum
     &           ,isumlo0,isumlo1
     &           ,isumhi0,isumhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer deriv_dir
      integer side
      REAL_T h(0:1)
      integer icoeflo0,icoeflo1
      integer icoefhi0,icoefhi1
      REAL_T coef(
     &           icoeflo0:icoefhi0,
     &           icoeflo1:icoefhi1)
      integer global(0:1)
      integer isumlo0,isumlo1
      integer isumhi0,isumhi1
      REAL_T sum(
     &           isumlo0:isumhi0,
     &           isumlo1:isumhi1)
      integer i0,i1, ii0,ii1, it0,it1, ic0,ic1, lo0,lo1,
     &     hi0,hi1, nlo0,nlo1, nhi0,nhi1, tii0,tii1,
     &     tdir, m, n, n_start, n_stop, t_start, t_stop
#if CH_SPACEDIM==3
     &     , tdir_other, t_start_other, t_stop_other
#endif
      REAL_T hfac, n_stencil(0:3), t_stencil_1(0:CH_SPACEDIM-1,0:4), t_stencil_2(0:CH_SPACEDIM-1,0:4),
     &     d, trans_grad_d, t_stencil_1_prod, t_stencil_2_prod
      
      ic0 = (isumlo0+isumhi0)/2
      ic1 = (isumlo1+isumhi1)/2
      hfac = dble(1 - 2*side) / h(deriv_dir)
      do tdir = 0, CH_SPACEDIM-1
         if (tdir .ne. dir) then
            hfac = hfac * h(tdir)
         endif
      enddo
      
      ii0=CHF_ID(0,dir)

      ii1=CHF_ID(1,dir)

      d = coef(global(0)+side*ii0,global(1)+side*ii1)
      if (dir .eq. deriv_dir) then
         n_stencil(0) =    one / 24.d0
         n_stencil(1) = -27.d0 / 24.d0
         n_stencil(2) =  27.d0 / 24.d0
         n_stencil(3) =   -one / 24.d0
      else
         n_stencil(0) = -one / 16.d0
         n_stencil(1) = nine / 16.d0
         n_stencil(2) = nine / 16.d0
         n_stencil(3) = -one / 16.d0
      endif
      n_start = side - 2
      n_stop = n_start + 3
      
      nlo0 = ic0 + ii0*n_start
      nhi0 = ic0 + ii0*n_stop
      nlo1 = ic1 + ii1*n_start
      nhi1 = ic1 + ii1*n_stop
      do m = 0, 4
         do n = 0, CH_SPACEDIM-1
            t_stencil_1(n,m) = zero
            t_stencil_2(n,m) = zero
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
         do tdir = 0, CH_SPACEDIM-1
            if (tdir .ne. dir) then
               do n = 0, CH_SPACEDIM-1
                  t_stencil_2(n,0) =  zero
               enddo
               t_stencil_2(tdir,0) =  one
               t_stencil_2(tdir,1) = -two
               t_stencil_2(tdir,2) =  one
               
      tii0=CHF_ID(0,tdir)

      tii1=CHF_ID(1,tdir)

               
               lo0 = nlo0 + tii0*t_start
               hi0 = nhi0 + tii0*t_stop
               lo1 = nlo1 + tii1*t_start
               hi1 = nhi1 + tii1*t_stop
               
                  do i1 = lo1, hi1
               do i0 = lo0, hi0
                        n = ii0*(i0-lo0) + ii1*(i1-lo1)
                        t_stencil_2_prod = 0 + ii0*t_stencil_2(1,i1-lo1) + ii1*t_stencil_2(0,i0-lo0)
                        sum(i0,i1) = sum(i0,i1)
     &                       + hfac * d * n_stencil(n) * t_stencil_2_prod / 24.d0
                        
                  enddo
                     enddo
            endif
         enddo
      else
         do n = 0, CH_SPACEDIM-1
            t_stencil_1(n,0) = zero
            t_stencil_2(n,0) = zero
         enddo
         t_stencil_1(deriv_dir,0) =    one / twelve
         t_stencil_1(deriv_dir,1) = -eight / twelve
         t_stencil_1(deriv_dir,2) =   zero
         t_stencil_1(deriv_dir,3) =  eight / twelve
         t_stencil_1(deriv_dir,4) =   -one / twelve
         t_stencil_2(deriv_dir,0) = -half
         t_stencil_2(deriv_dir,1) =   one
         t_stencil_2(deriv_dir,2) =  zero
         t_stencil_2(deriv_dir,3) =  -one
         t_stencil_2(deriv_dir,4) =  half
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
                  t_stencil_1_prod = 0 + ii0*t_stencil_1(1,i1-lo1) + ii1*t_stencil_1(0,i0-lo0)
                  t_stencil_2_prod = 0 + ii0*t_stencil_2(1,i1-lo1) + ii1*t_stencil_2(0,i0-lo0)
                  sum(i0,i1) = sum(i0,i1)
     &                 + hfac * d * n_stencil(n) * (t_stencil_1_prod + t_stencil_2_prod / 24.d0)
                  
            enddo
               enddo
#if CH_SPACEDIM==3
         do n = 0, CH_SPACEDIM-1
            t_stencil_2(n,0) = one
         enddo
         do n = 0, CH_SPACEDIM-1
            if (n .ne. dir .and. n .ne. deriv_dir) then
               tdir_other = n
            endif
         enddo
         t_stencil_2(deriv_dir,0) =    one / twelve
         t_stencil_2(deriv_dir,1) = -eight / twelve
         t_stencil_2(deriv_dir,2) =   zero
         t_stencil_2(deriv_dir,3) =  eight / twelve
         t_stencil_2(deriv_dir,4) =   -one / twelve
         t_stencil_2(tdir_other,0) =  one
         t_stencil_2(tdir_other,1) = -two
         t_stencil_2(tdir_other,2) =  one
         t_start_other = -1
         t_stop_other = 1
         
      tii0=CHF_ID(0,tdir_other)

      tii1=CHF_ID(1,tdir_other)

         
         lo0 = lo0 + tii0*t_start_other
         hi0 = hi0 + tii0*t_stop_other
         lo1 = lo1 + tii1*t_start_other
         hi1 = hi1 + tii1*t_stop_other
         
            do i1 = lo1, hi1
         do i0 = lo0, hi0
                  n = ii0*(i0-lo0) + ii1*(i1-lo1)
                  t_stencil_2_prod = ii0*t_stencil_2(1,i1-lo1)*t_stencil_2(2,i2-lo2)
     &                             + ii1*t_stencil_2(0,i0-lo0)*t_stencil_2(2,i2-lo2)
     &                             + ii2*t_stencil_2(0,i0-lo0)*t_stencil_2(1,i1-lo1)
                  sum(i0,i1) = sum(i0,i1)
     &                 + hfac * d * n_stencil(n) * t_stencil_2_prod / 24.d0
                  
            enddo
               enddo
#endif
      endif
      if (dir .eq. deriv_dir) then
         n_stencil(0) = -one
         n_stencil(1) =  one
      else
         n_stencil(0) = half
         n_stencil(1) = half
      endif
      n_start = side - 1
      n_stop = n_start + 1
      
      nlo0 = ic0 + ii0*n_start
      nhi0 = ic0 + ii0*n_stop
      nlo1 = ic1 + ii1*n_start
      nhi1 = ic1 + ii1*n_stop
      t_start = -1
      t_stop = 1
      do tdir = 0, CH_SPACEDIM-1
         if (tdir .ne. dir) then
            do n = 0, CH_SPACEDIM-1
               t_stencil_2(n,0) = zero
            enddo
            t_stencil_2(tdir,0) = -half
            t_stencil_2(tdir,1) =  zero
            t_stencil_2(tdir,2) =  half
            
      tii0=CHF_ID(0,tdir)

      tii1=CHF_ID(1,tdir)

            
            lo0 = nlo0 + tii0*t_start
            hi0 = nhi0 + tii0*t_stop
            lo1 = nlo1 + tii1*t_start
            hi1 = nhi1 + tii1*t_stop
            
      it0=CHF_ID(0,tdir)

      it1=CHF_ID(1,tdir)

            trans_grad_d = half * (
     &           coef(global(0)+ii0*side+it0,global(1)+ii1*side+it1)
     &        -  coef(global(0)+ii0*side-it0,global(1)+ii1*side-it1)
     &                            )
            
               do i1 = lo1, hi1
            do i0 = lo0, hi0
                     n = ii0*(i0-lo0) + ii1*(i1-lo1)
                     t_stencil_2_prod = 0 + ii0*t_stencil_2(1,i1-lo1) + ii1*t_stencil_2(0,i0-lo0)
                     sum(i0,i1) = sum(i0,i1)
     &                    + hfac * trans_grad_d * n_stencil(n) * t_stencil_2_prod / twelve
                     
               enddo
                  enddo
         endif
      enddo
      return
      end
      subroutine ACCUM_FLUX_STENCIL2_2D(
     &           dir
     &           ,deriv_dir
     &           ,side
     &           ,h
     &           ,coef
     &           ,icoeflo0,icoeflo1
     &           ,icoefhi0,icoefhi1
     &           ,global
     &           ,sum
     &           ,isumlo0,isumlo1
     &           ,isumhi0,isumhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer deriv_dir
      integer side
      REAL_T h(0:1)
      integer icoeflo0,icoeflo1
      integer icoefhi0,icoefhi1
      REAL_T coef(
     &           icoeflo0:icoefhi0,
     &           icoeflo1:icoefhi1)
      integer global(0:1)
      integer isumlo0,isumlo1
      integer isumhi0,isumhi1
      REAL_T sum(
     &           isumlo0:isumhi0,
     &           isumlo1:isumhi1)
      integer i0,i1, ii0,ii1, ic0,ic1, tii0,tii1,
     &   lo0,lo1, hi0,hi1, nlo0,nlo1, nhi0,nhi1,
     &   tdir, m, n, n_start, n_stop, t_start, t_stop
      REAL_T hfac, d, n_stencil(0:1), t_stencil(0:CH_SPACEDIM-1,0:2), t_stencil_prod
      
      ic0 = (isumlo0+isumhi0)/2
      ic1 = (isumlo1+isumhi1)/2
      hfac = (1 - 2*side)/ h(deriv_dir)
      do tdir = 0, CH_SPACEDIM-1
         if (tdir .ne. dir) then
            hfac = hfac * h(tdir)
         endif
      enddo
      
      ii0=CHF_ID(0,dir)

      ii1=CHF_ID(1,dir)

      d = coef(global(0)+side*ii0,global(1)+side*ii1)
      if (dir .eq. deriv_dir) then
         n_stencil(0) = -one
         n_stencil(1) =  one
      else
         n_stencil(0) = half
         n_stencil(1) = half
      endif
      n_start = side - 1
      n_stop = n_start + 1
      
      nlo0 = ic0 + ii0*n_start
      nhi0 = ic0 + ii0*n_stop
      nlo1 = ic1 + ii1*n_start
      nhi1 = ic1 + ii1*n_stop
      do m = 0, 2
         do n = 0, CH_SPACEDIM-1
            t_stencil(n,m) = zero
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
         do n = 0, CH_SPACEDIM-1
            t_stencil(n,0) = zero
         enddo
         t_stencil(deriv_dir,0) = -half
         t_stencil(deriv_dir,1) =  zero
         t_stencil(deriv_dir,2) =  half
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
                  t_stencil_prod = 0 + ii0*t_stencil(1,i1-lo1) + ii1*t_stencil(0,i0-lo0)
                  sum(i0,i1) = sum(i0,i1) + hfac * d * n_stencil(n) * t_stencil_prod
                  
            enddo
               enddo
      endif
      return
      end
