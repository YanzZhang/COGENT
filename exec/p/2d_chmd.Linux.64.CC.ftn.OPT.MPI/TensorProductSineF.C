#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SET_TENSORPRODUCTSINE_2D(
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
     &           ,t
     &           ,amp
     &           ,mode
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
      REAL_T t
      REAL_T amp
      REAL_T mode(0:1)
      integer i,j
      integer d,n
      REAL_T product, x, dsin
      do n=0,nfcomp-1
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        product = one
        do d=0,ncoordcomp-1
          x = coord(i,j,d)
          if (d.eq.0) then
            x = x - t
          endif
          product = product * dcos( two * Pi * mode(d) * x )
        enddo
        f(i,j,n) = amp * product
      
      enddo
      enddo
      enddo
      return
      end
