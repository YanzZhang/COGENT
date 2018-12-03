#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine REDISTRIBUTENEGATIVES_4D(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2,iphilo3
     &           ,iphihi0,iphihi1,iphihi2,iphihi3
     &           ,nphicomp
     &           ,iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
     &           ,iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
     &           ,ineighborBoxlo0,ineighborBoxlo1,ineighborBoxlo2,ineighborBoxlo3
     &           ,ineighborBoxhi0,ineighborBoxhi1,ineighborBoxhi2,ineighborBoxhi3
     &           ,refVal
     &           ,unable
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           iphilo3:iphihi3,
     &           0:nphicomp-1)
      integer iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
      integer iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
      integer ineighborBoxlo0,ineighborBoxlo1,ineighborBoxlo2,ineighborBoxlo3
      integer ineighborBoxhi0,ineighborBoxhi1,ineighborBoxhi2,ineighborBoxhi3
      REAL_T refVal
      integer unable
      integer n
      integer i,j,k,l
      integer m1,m2,m3,m4
      REAL_T deltaPhi
      REAL_T xisum,scale
      REAL_T xi(-3:3,-3:3,-3:3,-3:3)
      do n=0, (nphicomp-1)
        
      do l = iinteriorBoxlo3,iinteriorBoxhi3
      do k = iinteriorBoxlo2,iinteriorBoxhi2
      do j = iinteriorBoxlo1,iinteriorBoxhi1
      do i = iinteriorBoxlo0,iinteriorBoxhi0

          if (phi(i,j,k,l,n).lt.zero) then
            deltaPhi = -phi(i,j,k,l,n)
            if ((refVal+deltaPhi).ne.refVal) then
              xisum = zero
              
      do m4 = ineighborBoxlo3,ineighborBoxhi3
      do m3 = ineighborBoxlo2,ineighborBoxhi2
      do m2 = ineighborBoxlo1,ineighborBoxhi1
      do m1 = ineighborBoxlo0,ineighborBoxhi0

                 xi(m1,m2,m3,m4)
     &             = max(zero,phi(i+m1,j+m2,k+m3,l+m4,n))
                 xisum = xisum + xi(m1,m2,m3,m4)
              
      enddo
      enddo
      enddo
      enddo
              if (xisum.ge.deltaPhi) then
                scale = deltaPhi / xisum
                
      do m4 = ineighborBoxlo3,ineighborBoxhi3
      do m3 = ineighborBoxlo2,ineighborBoxhi2
      do m2 = ineighborBoxlo1,ineighborBoxhi1
      do m1 = ineighborBoxlo0,ineighborBoxhi0

                   xi(m1,m2,m3,m4)
     &               = xi(m1,m2,m3,m4) * scale
                
      enddo
      enddo
      enddo
      enddo
                
      do m4 = ineighborBoxlo3,ineighborBoxhi3
      do m3 = ineighborBoxlo2,ineighborBoxhi2
      do m2 = ineighborBoxlo1,ineighborBoxhi1
      do m1 = ineighborBoxlo0,ineighborBoxhi0

                   phi(i+m1,j+m2,k+m3,l+m4,n)
     &               = phi(i+m1,j+m2,k+m3,l+m4,n)
     &                   - xi(m1,m2,m3,m4)
                
      enddo
      enddo
      enddo
      enddo
                phi(i,j,k,l,n) = zero
              else
                unable = unable + 1
              endif
            else
              phi(i,j,k,l,n) = zero
            endif
          endif
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FINDANYNEGATIVES_4D(
     &           count
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2,iphilo3
     &           ,iphihi0,iphihi1,iphihi2,iphihi3
     &           ,nphicomp
     &           ,iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
     &           ,iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer count
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           iphilo3:iphihi3,
     &           0:nphicomp-1)
      integer iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
      integer iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
      integer n
      integer i,j,k,l
      REAL_T val
      do n=0, (nphicomp-1)
        
      do l = iinteriorBoxlo3,iinteriorBoxhi3
      do k = iinteriorBoxlo2,iinteriorBoxhi2
      do j = iinteriorBoxlo1,iinteriorBoxhi1
      do i = iinteriorBoxlo0,iinteriorBoxhi0

          if (phi(i,j,k,l,n).lt.zero) then
            count = count + 1
            return
          endif
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FINDALLNEGATIVES_4D(
     &           count
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2,iphilo3
     &           ,iphihi0,iphihi1,iphihi2,iphihi3
     &           ,nphicomp
     &           ,iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
     &           ,iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer count
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           iphilo3:iphihi3,
     &           0:nphicomp-1)
      integer iinteriorBoxlo0,iinteriorBoxlo1,iinteriorBoxlo2,iinteriorBoxlo3
      integer iinteriorBoxhi0,iinteriorBoxhi1,iinteriorBoxhi2,iinteriorBoxhi3
      integer n
      integer i,j,k,l
      do n=0, (nphicomp-1)
        
      do l = iinteriorBoxlo3,iinteriorBoxhi3
      do k = iinteriorBoxlo2,iinteriorBoxhi2
      do j = iinteriorBoxlo1,iinteriorBoxhi1
      do i = iinteriorBoxlo0,iinteriorBoxhi0

          if (phi(i,j,k,l,n).lt.zero) then
            count = count + 1
          endif
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
