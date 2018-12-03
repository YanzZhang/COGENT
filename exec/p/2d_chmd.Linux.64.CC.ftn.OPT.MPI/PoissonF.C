#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine COMPUTE_LAPLACIAN_COEFFICIENTS_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,n
     &           ,inlo0,inlo1
     &           ,inhi0,inhi1
     &           ,nncomp
     &           ,njinverse
     &           ,injinverselo0,injinverselo1
     &           ,injinversehi0,injinversehi1
     &           ,nnjinversecomp
     &           ,coef
     &           ,icoeflo0,icoeflo1
     &           ,icoefhi0,icoefhi1
     &           ,ncoefcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer nncomp
      integer inlo0,inlo1
      integer inhi0,inhi1
      REAL_T n(
     &           inlo0:inhi0,
     &           inlo1:inhi1,
     &           0:nncomp-1)
      integer nnjinversecomp
      integer injinverselo0,injinverselo1
      integer injinversehi0,injinversehi1
      REAL_T njinverse(
     &           injinverselo0:injinversehi0,
     &           injinverselo1:injinversehi1,
     &           0:nnjinversecomp-1)
      integer ncoefcomp
      integer icoeflo0,icoeflo1
      integer icoefhi0,icoefhi1
      REAL_T coef(
     &           icoeflo0:icoefhi0,
     &           icoeflo1:icoefhi1,
     &           0:ncoefcomp-1)
      integer i,j, row, col, m
      double precision n_mat(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1),
     &       nji_mat(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1),
     &       d_mat(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1),
     &       dnji_mat(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1),
     &       coef_mat(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==3
         d_mat(0,0) = one
         d_mat(0,1) = zero
         d_mat(0,2) = zero
         d_mat(1,0) = zero
         d_mat(1,1) = one
         d_mat(1,2) = zero
         d_mat(2,0) = zero
         d_mat(2,1) = zero
         d_mat(2,2) = one
         nji_mat(0,0) = njinverse(i,j,0)
         nji_mat(0,1) = njinverse(i,j,1)
         nji_mat(0,2) = njinverse(i,j,2)
         nji_mat(1,0) = njinverse(i,j,3)
         nji_mat(1,1) = njinverse(i,j,4)
         nji_mat(1,2) = njinverse(i,j,5)
         nji_mat(2,0) = njinverse(i,j,6)
         nji_mat(2,1) = njinverse(i,j,7)
         nji_mat(2,2) = njinverse(i,j,8)
         n_mat(0,0) = n(i,j,0)
         n_mat(0,1) = n(i,j,1)
         n_mat(0,2) = n(i,j,2)
         n_mat(1,0) = n(i,j,3)
         n_mat(1,1) = n(i,j,4)
         n_mat(1,2) = n(i,j,5)
         n_mat(2,0) = n(i,j,6)
         n_mat(2,1) = n(i,j,7)
         n_mat(2,2) = n(i,j,8)
#else
         d_mat(0,0) = one
         d_mat(0,1) = zero
         d_mat(1,0) = zero
         d_mat(1,1) = one
         nji_mat(0,0) = njinverse(i,j,0)
         nji_mat(0,1) = njinverse(i,j,1)
         nji_mat(1,0) = njinverse(i,j,2)
         nji_mat(1,1) = njinverse(i,j,3)
         n_mat(0,0) = n(i,j,0)
         n_mat(0,1) = n(i,j,1)
         n_mat(1,0) = n(i,j,2)
         n_mat(1,1) = n(i,j,3)
#endif
         do row = 0, CH_SPACEDIM-1
            do col = 0, CH_SPACEDIM-1
               dnji_mat(row,col) = zero
               do m = 0, CH_SPACEDIM-1
                  dnji_mat(row,col) = dnji_mat(row,col) + d_mat(row,m) * nji_mat(m,col)
               enddo
            enddo
         enddo
         do row = 0, CH_SPACEDIM-1
            do col = 0, CH_SPACEDIM-1
               coef_mat(row,col) = zero
               do m = 0, CH_SPACEDIM-1
                  coef_mat(row,col) = coef_mat(row,col) + n_mat(m,row) * dnji_mat(m,col)
               enddo
            enddo
         enddo
#if CH_SPACEDIM==3
         coef(i,j,0) = coef_mat(0,0)
         coef(i,j,1) = coef_mat(0,1)
         coef(i,j,2) = coef_mat(0,2)
         coef(i,j,3) = coef_mat(1,0)
         coef(i,j,4) = coef_mat(1,1)
         coef(i,j,5) = coef_mat(1,2)
         coef(i,j,6) = coef_mat(2,0)
         coef(i,j,7) = coef_mat(2,1)
         coef(i,j,8) = coef_mat(2,2)
#else
         coef(i,j,0) = coef_mat(0,0)
         coef(i,j,1) = coef_mat(0,1)
         coef(i,j,2) = coef_mat(1,0)
         coef(i,j,3) = coef_mat(1,1)
#endif
      
      enddo
      enddo
      return
      end
