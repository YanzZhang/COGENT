#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine COMPUTE_MAPPED_COEFFICIENTS_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,unmapped_coef
     &           ,iunmapped_coeflo0,iunmapped_coeflo1
     &           ,iunmapped_coefhi0,iunmapped_coefhi1
     &           ,nunmapped_coefcomp
     &           ,n
     &           ,inlo0,inlo1
     &           ,inhi0,inhi1
     &           ,nncomp
     &           ,jinverse
     &           ,ijinverselo0,ijinverselo1
     &           ,ijinversehi0,ijinversehi1
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
      integer nunmapped_coefcomp
      integer iunmapped_coeflo0,iunmapped_coeflo1
      integer iunmapped_coefhi0,iunmapped_coefhi1
      REAL_T unmapped_coef(
     &           iunmapped_coeflo0:iunmapped_coefhi0,
     &           iunmapped_coeflo1:iunmapped_coefhi1,
     &           0:nunmapped_coefcomp-1)
      integer nncomp
      integer inlo0,inlo1
      integer inhi0,inhi1
      REAL_T n(
     &           inlo0:inhi0,
     &           inlo1:inhi1,
     &           0:nncomp-1)
      integer ijinverselo0,ijinverselo1
      integer ijinversehi0,ijinversehi1
      REAL_T jinverse(
     &           ijinverselo0:ijinversehi0,
     &           ijinverselo1:ijinversehi1)
      integer ncoefcomp
      integer icoeflo0,icoeflo1
      integer icoefhi0,icoefhi1
      REAL_T coef(
     &           icoeflo0:icoefhi0,
     &           icoeflo1:icoefhi1,
     &           0:ncoefcomp-1)
      integer i,j, row, col, m
      double precision n_mat(CH_SPACEDIM,CH_SPACEDIM), nji_mat(CH_SPACEDIM,CH_SPACEDIM),
     &       d_mat(CH_SPACEDIM,CH_SPACEDIM), dnji_mat(CH_SPACEDIM,CH_SPACEDIM),
     &       coef_mat(CH_SPACEDIM,CH_SPACEDIM)
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

         m = -1
         do row = 1, CH_SPACEDIM
            do col = 1, CH_SPACEDIM
               m = m + 1
               d_mat(row,col) = unmapped_coef(i,j,m)
               nji_mat(row,col) = n(i,j,m) * jinverse(i,j)
            enddo
         enddo
         do row = 1, CH_SPACEDIM
            do col = 1, CH_SPACEDIM
               dnji_mat(row,col) = zero
               do m = 1, CH_SPACEDIM
                  dnji_mat(row,col) = dnji_mat(row,col) + d_mat(row,m) * nji_mat(m,col)
               enddo
            enddo
         enddo
         m = -1
         do row = 1, CH_SPACEDIM
            do col = 1, CH_SPACEDIM
               m = m + 1
               n_mat(row,col) = n(i,j,m)
            enddo
         enddo
         do row = 1, CH_SPACEDIM
            do col = 1, CH_SPACEDIM
               coef_mat(row,col) = zero
               do m = 1, CH_SPACEDIM
                  coef_mat(row,col) = coef_mat(row,col) + n_mat(m,row) * dnji_mat(m,col)
               enddo
            enddo
         enddo
         m = -1
         do row = 1, CH_SPACEDIM
            do col = 1, CH_SPACEDIM
               m = m + 1
               coef(i,j,m) = coef_mat(row,col)
            enddo
         enddo
      
      enddo
      enddo
      return
      end
