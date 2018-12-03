#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine COMPUTE_ELECTRON_DENSITY_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,te
     &           ,itelo0,itelo1
     &           ,itehi0,itehi1
     &           ,ne
     &           ,inelo0,inelo1
     &           ,inehi0,inehi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      integer itelo0,itelo1
      integer itehi0,itehi1
      REAL_T te(
     &           itelo0:itehi0,
     &           itelo1:itehi1)
      integer inelo0,inelo1
      integer inehi0,inehi1
      REAL_T ne(
     &           inelo0:inehi0,
     &           inelo1:inehi1)
      integer i,j
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

         ne(i,j) = dexp(phi(i,j)/te(i,j))
      
      enddo
      enddo
      return
      end
