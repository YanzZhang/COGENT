#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine INSPECT_3D(
     &           igridboxlo0,igridboxlo1,igridboxlo2
     &           ,igridboxhi0,igridboxhi1,igridboxhi2
     &           ,data
     &           ,idatalo0,idatalo1,idatalo2
     &           ,idatahi0,idatahi1,idatahi2
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2
      integer igridboxhi0,igridboxhi1,igridboxhi2
      integer ndatacomp
      integer idatalo0,idatalo1,idatalo2
      integer idatahi0,idatahi1,idatahi2
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           idatalo2:idatahi2,
     &           0:ndatacomp-1)
      integer i,j,k,n
      do n=0,ndatacomp-1
      
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

      data(i,j,k,n)=data(i,j,k,n)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
