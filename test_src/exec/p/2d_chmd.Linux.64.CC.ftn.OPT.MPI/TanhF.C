#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SET_TANH_2D(
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
     &           ,f1
     &           ,x1
     &           ,f2
     &           ,x2
     &           ,x0
     &           ,width
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
      REAL_T f1
      REAL_T x1
      REAL_T f2
      REAL_T x2
      REAL_T x0
      REAL_T width
      integer i,j
      integer n
      REAL_T a,b,arg,invwidth
      REAL_T tanhx1,tanhx2
      REAL_T tanh
      invwidth = one / width
      tanhx1 = tanh( (x1 - x0) * invwidth )
      tanhx2 = tanh( (x2 - x0) * invwidth )
      a = (f1 - f2) / (tanhx1 - tanhx2)
      b = f2 - a * tanhx2
      do n=0,nphicomp-1
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        x(i,j,0) = (i + half)*dx(0)
        arg = (x(i,j,0) - x0) * invwidth
        phi(i,j,n) = a * tanh( arg ) + b
      
      enddo
      enddo
      enddo
      return
      end
