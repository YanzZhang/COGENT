#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SET_LOCALIZED_2D(
     &           f
     &           ,iflo0,iflo1
     &           ,ifhi0,ifhi1
     &           ,nfcomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,coord
     &           ,icoordlo0,icoordlo1
     &           ,icoordhi0,icoordhi1
     &           ,ncoordcomp
     &           ,amp
     &           ,center
     &           ,width
     &           ,floor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfcomp
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           0:nfcomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ncoordcomp
      integer icoordlo0,icoordlo1
      integer icoordhi0,icoordhi1
      REAL_T coord(
     &           icoordlo0:icoordhi0,
     &           icoordlo1:icoordhi1,
     &           0:ncoordcomp-1)
      REAL_T amp
      REAL_T center(0:1)
      REAL_T width(0:1)
      REAL_T floor
      integer i,j
      integer d,n
      REAL_T beta, x
      do n=0,nfcomp-1
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        beta = zero
        do d=0,ncoordcomp-1
          x = coord(i,j,d) - center(d)
          beta = beta + ( x / width(d) )**2
        enddo
        f(i,j,n) = amp * dexp(-beta) + floor
      
      enddo
      enddo
      enddo
      return
      end
