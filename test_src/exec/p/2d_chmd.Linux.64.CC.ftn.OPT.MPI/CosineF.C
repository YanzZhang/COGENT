#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SET_COSINE_2D(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,dx
     &           ,shift
     &           ,constant
     &           ,amplitude
     &           ,mode
     &           ,phase
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL_T x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      REAL_T dx(0:1)
      REAL_T shift
      REAL_T constant
      REAL_T amplitude
      REAL_T mode(0:1)
      REAL_T phase(0:1)
      integer i,j
      integer d,n
      REAL_T arg,product
      do n=0,nphicomp-1
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        x(i,j,0) = (i + half)*dx(0) + shift
        product = one
        do d=0,nxcomp-1
          arg = mode(d) * x(i,j,d) + phase(d)
          product = product * cos( arg )
        enddo
        phi(i,j,n) = constant + amplitude * product
      
      enddo
      enddo
      enddo
      return
      end
