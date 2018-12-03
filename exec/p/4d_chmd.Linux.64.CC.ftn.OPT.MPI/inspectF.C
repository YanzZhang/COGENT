#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine INSPECT_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,data
     &           ,idatalo0,idatalo1,idatalo2,idatalo3
     &           ,idatahi0,idatahi1,idatahi2,idatahi3
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer ndatacomp
      integer idatalo0,idatalo1,idatalo2,idatalo3
      integer idatahi0,idatahi1,idatahi2,idatahi3
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           idatalo2:idatahi2,
     &           idatalo3:idatahi3,
     &           0:ndatacomp-1)
      integer i,j,k,l,n
      do n=0,ndatacomp-1
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

      data(i,j,k,l,n)=data(i,j,k,l,n)
      
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
