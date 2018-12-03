#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine INSPECT_2D(
     &           igridboxlo0,igridboxlo1
     &           ,igridboxhi0,igridboxhi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           0:ndatacomp-1)
      integer i,j,n
      do n=0,ndatacomp-1
      
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

      data(i,j,n)=data(i,j,n)
      
      enddo
      enddo
      enddo
      return
      end
