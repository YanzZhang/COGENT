#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine CENTEREDLAPLACIAN_4D(
     &           lapPhi
     &           ,ilapPhilo0,ilapPhilo1,ilapPhilo2,ilapPhilo3
     &           ,ilapPhihi0,ilapPhihi1,ilapPhihi2,ilapPhihi3
     &           ,nlapPhicomp
     &           ,facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
     &           ,ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
     &           ,idir
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlapPhicomp
      integer ilapPhilo0,ilapPhilo1,ilapPhilo2,ilapPhilo3
      integer ilapPhihi0,ilapPhihi1,ilapPhihi2,ilapPhihi3
      REAL_T lapPhi(
     &           ilapPhilo0:ilapPhihi0,
     &           ilapPhilo1:ilapPhihi1,
     &           ilapPhilo2:ilapPhihi2,
     &           ilapPhilo3:ilapPhihi3,
     &           0:nlapPhicomp-1)
      integer nfacePhicomp
      integer ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
      integer ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
      REAL_T facePhi(
     &           ifacePhilo0:ifacePhihi0,
     &           ifacePhilo1:ifacePhihi1,
     &           ifacePhilo2:ifacePhihi2,
     &           ifacePhilo3:ifacePhihi3,
     &           0:nfacePhicomp-1)
      integer ncellPhicomp
      integer icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
      integer icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
      REAL_T cellPhi(
     &           icellPhilo0:icellPhihi0,
     &           icellPhilo1:icellPhihi1,
     &           icellPhilo2:icellPhihi2,
     &           icellPhilo3:icellPhihi3,
     &           0:ncellPhicomp-1)
      integer ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
      integer ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
      integer idir
      REAL_T dx
      integer n, i,j,k,l
      integer ii,jj,kk,ll
      REAL_T factor
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      factor = three/(dx*dx)
      do n=0, (nfacePhicomp-1)
         
      do l = ifaceBoxlo3,ifaceBoxhi3
      do k = ifaceBoxlo2,ifaceBoxhi2
      do j = ifaceBoxlo1,ifaceBoxhi1
      do i = ifaceBoxlo0,ifaceBoxhi0

            lapPhi(i,j,k,l,n)
     &        = lapPhi(i,j,k,l,n)
     &        + factor * ( cellPhi(i,j,k,l,n)
     &                   + cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &                   - two * facePhi(i,j,k,l,n) )
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine CCLAPLACIAN_4D(
     &           lapPhi
     &           ,ilapPhilo0,ilapPhilo1,ilapPhilo2,ilapPhilo3
     &           ,ilapPhihi0,ilapPhihi1,ilapPhihi2,ilapPhihi3
     &           ,nlapPhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,icellBoxlo0,icellBoxlo1,icellBoxlo2,icellBoxlo3
     &           ,icellBoxhi0,icellBoxhi1,icellBoxhi2,icellBoxhi3
     &           ,idir
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlapPhicomp
      integer ilapPhilo0,ilapPhilo1,ilapPhilo2,ilapPhilo3
      integer ilapPhihi0,ilapPhihi1,ilapPhihi2,ilapPhihi3
      REAL_T lapPhi(
     &           ilapPhilo0:ilapPhihi0,
     &           ilapPhilo1:ilapPhihi1,
     &           ilapPhilo2:ilapPhihi2,
     &           ilapPhilo3:ilapPhihi3,
     &           0:nlapPhicomp-1)
      integer ncellPhicomp
      integer icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
      integer icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
      REAL_T cellPhi(
     &           icellPhilo0:icellPhihi0,
     &           icellPhilo1:icellPhihi1,
     &           icellPhilo2:icellPhihi2,
     &           icellPhilo3:icellPhihi3,
     &           0:ncellPhicomp-1)
      integer icellBoxlo0,icellBoxlo1,icellBoxlo2,icellBoxlo3
      integer icellBoxhi0,icellBoxhi1,icellBoxhi2,icellBoxhi3
      integer idir
      REAL_T dx
      integer n, i,j,k,l
      integer ii,jj,kk,ll
      REAL_T factor
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      factor = one/(dx*dx)
      do n=0, (nlapPhicomp-1)
         
      do l = icellBoxlo3,icellBoxhi3
      do k = icellBoxlo2,icellBoxhi2
      do j = icellBoxlo1,icellBoxhi1
      do i = icellBoxlo0,icellBoxhi0

            lapPhi(i,j,k,l,n)
     &        = lapPhi(i,j,k,l,n)
     &        + factor * ( cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &                   + cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &                   - two * cellPhi(i,j,k,l,n) )
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine LIMITFACEVALUES_4D(
     &           facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,centeredLap
     &           ,icenteredLaplo0,icenteredLaplo1,icenteredLaplo2,icenteredLaplo3
     &           ,icenteredLaphi0,icenteredLaphi1,icenteredLaphi2,icenteredLaphi3
     &           ,ncenteredLapcomp
     &           ,LRLap
     &           ,iLRLaplo0,iLRLaplo1,iLRLaplo2,iLRLaplo3
     &           ,iLRLaphi0,iLRLaphi1,iLRLaphi2,iLRLaphi3
     &           ,nLRLapcomp
     &           ,ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
     &           ,ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
     &           ,idir
     &           ,dx
     &           ,limitC
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfacePhicomp
      integer ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
      integer ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
      REAL_T facePhi(
     &           ifacePhilo0:ifacePhihi0,
     &           ifacePhilo1:ifacePhihi1,
     &           ifacePhilo2:ifacePhihi2,
     &           ifacePhilo3:ifacePhihi3,
     &           0:nfacePhicomp-1)
      integer ncellPhicomp
      integer icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
      integer icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
      REAL_T cellPhi(
     &           icellPhilo0:icellPhihi0,
     &           icellPhilo1:icellPhihi1,
     &           icellPhilo2:icellPhihi2,
     &           icellPhilo3:icellPhihi3,
     &           0:ncellPhicomp-1)
      integer ncenteredLapcomp
      integer icenteredLaplo0,icenteredLaplo1,icenteredLaplo2,icenteredLaplo3
      integer icenteredLaphi0,icenteredLaphi1,icenteredLaphi2,icenteredLaphi3
      REAL_T centeredLap(
     &           icenteredLaplo0:icenteredLaphi0,
     &           icenteredLaplo1:icenteredLaphi1,
     &           icenteredLaplo2:icenteredLaphi2,
     &           icenteredLaplo3:icenteredLaphi3,
     &           0:ncenteredLapcomp-1)
      integer nLRLapcomp
      integer iLRLaplo0,iLRLaplo1,iLRLaplo2,iLRLaplo3
      integer iLRLaphi0,iLRLaphi1,iLRLaphi2,iLRLaphi3
      REAL_T LRLap(
     &           iLRLaplo0:iLRLaphi0,
     &           iLRLaplo1:iLRLaphi1,
     &           iLRLaplo2:iLRLaphi2,
     &           iLRLaplo3:iLRLaphi3,
     &           0:nLRLapcomp-1)
      integer ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
      integer ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
      integer idir
      REAL_T dx
      REAL_T limitC
      integer n, i,j,k,l
      integer ii,jj,kk,ll
      REAL_T d2phi, d2C, d2L, d2R
      REAL_T h2On3
      REAL_T diffl, diffr
      REAL_T s
      logical sameSign
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      do n=0, (nfacePhicomp-1)
         
      do l = ifaceBoxlo3,ifaceBoxhi3
      do k = ifaceBoxlo2,ifaceBoxhi2
      do j = ifaceBoxlo1,ifaceBoxhi1
      do i = ifaceBoxlo0,ifaceBoxhi0

           diffl = facePhi(i,j,k,l,n)
     &           - cellPhi(i-ii,j-jj,k-kk,l-ll,n)
           diffr = cellPhi(i,j,k,l,n)
     &           - facePhi(i,j,k,l,n)
           if (diffl*diffr.lt.zero) then
              d2L = LRLap(i-ii,j-jj,k-kk,l-ll,n)
              d2R = LRLap(i,j,k,l,n)
              d2C = centeredLap(i,j,k,l,n)
              if ((d2R - d2C)*(d2C-d2L) .lt.zero) then
                 s = sign(one,d2C)
                 d2phi = s*max(min(limitC*s*d2L,limitC*s*d2R,s*d2C),zero)
                 facePhi(i,j,k,l,n)
     &                = half * ( cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &                      + cellPhi(i,j,k,l,n) )
     &                - sixth * d2phi * dx**2
              endif
           endif
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTEA6_4D(
     &           a6
     &           ,ia6lo0,ia6lo1,ia6lo2,ia6lo3
     &           ,ia6hi0,ia6hi1,ia6hi2,ia6hi3
     &           ,na6comp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,iccboxlo0,iccboxlo1,iccboxlo2,iccboxlo3
     &           ,iccboxhi0,iccboxhi1,iccboxhi2,iccboxhi3
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer na6comp
      integer ia6lo0,ia6lo1,ia6lo2,ia6lo3
      integer ia6hi0,ia6hi1,ia6hi2,ia6hi3
      REAL_T a6(
     &           ia6lo0:ia6hi0,
     &           ia6lo1:ia6hi1,
     &           ia6lo2:ia6hi2,
     &           ia6lo3:ia6hi3,
     &           0:na6comp-1)
      integer ncellPhicomp
      integer icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
      integer icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
      REAL_T cellPhi(
     &           icellPhilo0:icellPhihi0,
     &           icellPhilo1:icellPhihi1,
     &           icellPhilo2:icellPhihi2,
     &           icellPhilo3:icellPhihi3,
     &           0:ncellPhicomp-1)
      integer nfacePhicomp
      integer ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
      integer ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
      REAL_T facePhi(
     &           ifacePhilo0:ifacePhihi0,
     &           ifacePhilo1:ifacePhihi1,
     &           ifacePhilo2:ifacePhihi2,
     &           ifacePhilo3:ifacePhihi3,
     &           0:nfacePhicomp-1)
      integer iccboxlo0,iccboxlo1,iccboxlo2,iccboxlo3
      integer iccboxhi0,iccboxhi1,iccboxhi2,iccboxhi3
      integer dir
      integer i,j,k,l, n
      integer ii,jj,kk,ll
      ii = CHF_ID(dir,0)
                jj = CHF_ID(dir,1)
                kk = CHF_ID(dir,2)
                ll = CHF_ID(dir,3)
      do n=0, na6comp-1
         
      do l = iccboxlo3,iccboxhi3
      do k = iccboxlo2,iccboxhi2
      do j = iccboxlo1,iccboxhi1
      do i = iccboxlo0,iccboxhi0

           a6(i,j,k,l,n)
     &       = six * cellPhi(i,j,k,l,n)
     &       - three * ( facePhi(i,j,k,l,n)
     &                 + facePhi(i+ii,j+jj,k+kk,l+ll,n) )
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine LIMITD2A_4D(
     &           D2PhiLim
     &           ,iD2PhiLimlo0,iD2PhiLimlo1,iD2PhiLimlo2,iD2PhiLimlo3
     &           ,iD2PhiLimhi0,iD2PhiLimhi1,iD2PhiLimhi2,iD2PhiLimhi3
     &           ,nD2PhiLimcomp
     &           ,D2Phi
     &           ,iD2Philo0,iD2Philo1,iD2Philo2,iD2Philo3
     &           ,iD2Phihi0,iD2Phihi1,iD2Phihi2,iD2Phihi3
     &           ,nD2Phicomp
     &           ,lapPhi
     &           ,ilapPhilo0,ilapPhilo1,ilapPhilo2,ilapPhilo3
     &           ,ilapPhihi0,ilapPhihi1,ilapPhihi2,ilapPhihi3
     &           ,nlapPhicomp
     &           ,icellBoxlo0,icellBoxlo1,icellBoxlo2,icellBoxlo3
     &           ,icellBoxhi0,icellBoxhi1,icellBoxhi2,icellBoxhi3
     &           ,dir
     &           ,limitC
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nD2PhiLimcomp
      integer iD2PhiLimlo0,iD2PhiLimlo1,iD2PhiLimlo2,iD2PhiLimlo3
      integer iD2PhiLimhi0,iD2PhiLimhi1,iD2PhiLimhi2,iD2PhiLimhi3
      REAL_T D2PhiLim(
     &           iD2PhiLimlo0:iD2PhiLimhi0,
     &           iD2PhiLimlo1:iD2PhiLimhi1,
     &           iD2PhiLimlo2:iD2PhiLimhi2,
     &           iD2PhiLimlo3:iD2PhiLimhi3,
     &           0:nD2PhiLimcomp-1)
      integer nD2Phicomp
      integer iD2Philo0,iD2Philo1,iD2Philo2,iD2Philo3
      integer iD2Phihi0,iD2Phihi1,iD2Phihi2,iD2Phihi3
      REAL_T D2Phi(
     &           iD2Philo0:iD2Phihi0,
     &           iD2Philo1:iD2Phihi1,
     &           iD2Philo2:iD2Phihi2,
     &           iD2Philo3:iD2Phihi3,
     &           0:nD2Phicomp-1)
      integer nlapPhicomp
      integer ilapPhilo0,ilapPhilo1,ilapPhilo2,ilapPhilo3
      integer ilapPhihi0,ilapPhihi1,ilapPhihi2,ilapPhihi3
      REAL_T lapPhi(
     &           ilapPhilo0:ilapPhihi0,
     &           ilapPhilo1:ilapPhihi1,
     &           ilapPhilo2:ilapPhihi2,
     &           ilapPhilo3:ilapPhihi3,
     &           0:nlapPhicomp-1)
      integer icellBoxlo0,icellBoxlo1,icellBoxlo2,icellBoxlo3
      integer icellBoxhi0,icellBoxhi1,icellBoxhi2,icellBoxhi3
      integer dir
      REAL_T limitC
      integer i,j,k,l,n
      integer ii,jj,kk,ll
      REAL_T d2al, d2ar, d2ac, d2a, minVal
      REAL_T s
      logical sameSign
      ii = CHF_ID(dir,0)
                jj = CHF_ID(dir,1)
                kk = CHF_ID(dir,2)
                ll = CHF_ID(dir,3)
      do n=0, nD2PhiLimcomp-1
         
      do l = icellBoxlo3,icellBoxhi3
      do k = icellBoxlo2,icellBoxhi2
      do j = icellBoxlo1,icellBoxhi1
      do i = icellBoxlo0,icellBoxhi0

           d2a  = D2Phi(i,j,k,l,n)
           d2al = lapPhi(i-ii,j-jj,k-kk,l-ll,n)
           d2ar = lapPhi(i+ii,j+jj,k+kk,l+ll,n)
           d2ac = lapPhi(i,j,k,l,n)
           if ((d2ac - d2al)*(d2ar - d2ac) .lt. zero) then
              s = sign(one,d2a)
              minVal = min(s*d2a,limitC*s*d2al,limitC*s*d2ar,limitC*s*d2ac)
              D2PhiLim(i,j,k,l,n) = max(minVal,zero)
           else
              D2PhiLim(i,j,k,l,n) = abs(d2a)
           endif
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine LEFTRIGHTSTATES_4D(
     &           phiLeft
     &           ,iphiLeftlo0,iphiLeftlo1,iphiLeftlo2,iphiLeftlo3
     &           ,iphiLefthi0,iphiLefthi1,iphiLefthi2,iphiLefthi3
     &           ,nphiLeftcomp
     &           ,phiRight
     &           ,iphiRightlo0,iphiRightlo1,iphiRightlo2,iphiRightlo3
     &           ,iphiRighthi0,iphiRighthi1,iphiRighthi2,iphiRighthi3
     &           ,nphiRightcomp
     &           ,facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,D2PhiLim
     &           ,iD2PhiLimlo0,iD2PhiLimlo1,iD2PhiLimlo2,iD2PhiLimlo3
     &           ,iD2PhiLimhi0,iD2PhiLimhi1,iD2PhiLimhi2,iD2PhiLimhi3
     &           ,nD2PhiLimcomp
     &           ,D2Phi
     &           ,iD2Philo0,iD2Philo1,iD2Philo2,iD2Philo3
     &           ,iD2Phihi0,iD2Phihi1,iD2Phihi2,iD2Phihi3
     &           ,nD2Phicomp
     &           ,icellBoxlo0,icellBoxlo1,icellBoxlo2,icellBoxlo3
     &           ,icellBoxhi0,icellBoxhi1,icellBoxhi2,icellBoxhi3
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphiLeftcomp
      integer iphiLeftlo0,iphiLeftlo1,iphiLeftlo2,iphiLeftlo3
      integer iphiLefthi0,iphiLefthi1,iphiLefthi2,iphiLefthi3
      REAL_T phiLeft(
     &           iphiLeftlo0:iphiLefthi0,
     &           iphiLeftlo1:iphiLefthi1,
     &           iphiLeftlo2:iphiLefthi2,
     &           iphiLeftlo3:iphiLefthi3,
     &           0:nphiLeftcomp-1)
      integer nphiRightcomp
      integer iphiRightlo0,iphiRightlo1,iphiRightlo2,iphiRightlo3
      integer iphiRighthi0,iphiRighthi1,iphiRighthi2,iphiRighthi3
      REAL_T phiRight(
     &           iphiRightlo0:iphiRighthi0,
     &           iphiRightlo1:iphiRighthi1,
     &           iphiRightlo2:iphiRighthi2,
     &           iphiRightlo3:iphiRighthi3,
     &           0:nphiRightcomp-1)
      integer nfacePhicomp
      integer ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
      integer ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
      REAL_T facePhi(
     &           ifacePhilo0:ifacePhihi0,
     &           ifacePhilo1:ifacePhihi1,
     &           ifacePhilo2:ifacePhihi2,
     &           ifacePhilo3:ifacePhihi3,
     &           0:nfacePhicomp-1)
      integer ncellPhicomp
      integer icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
      integer icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
      REAL_T cellPhi(
     &           icellPhilo0:icellPhihi0,
     &           icellPhilo1:icellPhihi1,
     &           icellPhilo2:icellPhihi2,
     &           icellPhilo3:icellPhihi3,
     &           0:ncellPhicomp-1)
      integer nD2PhiLimcomp
      integer iD2PhiLimlo0,iD2PhiLimlo1,iD2PhiLimlo2,iD2PhiLimlo3
      integer iD2PhiLimhi0,iD2PhiLimhi1,iD2PhiLimhi2,iD2PhiLimhi3
      REAL_T D2PhiLim(
     &           iD2PhiLimlo0:iD2PhiLimhi0,
     &           iD2PhiLimlo1:iD2PhiLimhi1,
     &           iD2PhiLimlo2:iD2PhiLimhi2,
     &           iD2PhiLimlo3:iD2PhiLimhi3,
     &           0:nD2PhiLimcomp-1)
      integer nD2Phicomp
      integer iD2Philo0,iD2Philo1,iD2Philo2,iD2Philo3
      integer iD2Phihi0,iD2Phihi1,iD2Phihi2,iD2Phihi3
      REAL_T D2Phi(
     &           iD2Philo0:iD2Phihi0,
     &           iD2Philo1:iD2Phihi1,
     &           iD2Philo2:iD2Phihi2,
     &           iD2Philo3:iD2Phihi3,
     &           0:nD2Phicomp-1)
      integer icellBoxlo0,icellBoxlo1,icellBoxlo2,icellBoxlo3
      integer icellBoxhi0,icellBoxhi1,icellBoxhi2,icellBoxhi3
      integer dir
      integer i,j,k,l,n
      integer ii,jj,kk,ll
      REAL_T ratio, leftVal, rightVal, cellVal
      REAL_T maxCheckFC, maxCheckCC
      REAL_T zeroTol
      REAL_T phil, phir
      REAL_T dphifacel, dphifacer, dphifacemin
      REAL_T dphibarl, dphibarr, dphibarmin
      REAL_T dphichkl, dphichkr
      REAL_T delphil, delphir
      REAL_T phimax, root, diff
      logical bigl,bigr,extremum
      zeroTol = 1.0d-10
      ii = CHF_ID(dir,0)
                jj = CHF_ID(dir,1)
                kk = CHF_ID(dir,2)
                ll = CHF_ID(dir,3)
      do n=0, nD2PhiLimcomp-1
        
      do l = icellBoxlo3,icellBoxhi3
      do k = icellBoxlo2,icellBoxhi2
      do j = icellBoxlo1,icellBoxhi1
      do i = icellBoxlo0,icellBoxhi0

          phil = facePhi(i,j,k,l,n)
     &         - cellPhi(i,j,k,l,n)
          phir = facePhi(i+ii,j+jj,k+kk,l+ll,n)
     &         - cellPhi(i,j,k,l,n)
          bigl = abs(phil).gt.2*abs(phir)
          bigr = abs(phir).gt.2*abs(phil)
          extremum = .false.
          if (phil*phir.ge.zero) then
            extremum = .true.
          else if (bigl.or.bigr) then
            dphifacel = facePhi(i,j,k,l,n)
     &                - facePhi(i-ii,j-jj,k-kk,l-ll,n)
            dphifacer = facePhi(i+2*ii,j+2*jj,k+2*kk,l+2*ll,n)
     &                - facePhi(i+ii,j+jj,k+kk,l+ll,n)
            dphifacemin = min(abs(dphifacel),abs(dphifacer))
            dphibarr = cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &               - cellPhi(i,j,k,l,n)
            dphibarl = cellPhi(i,j,k,l,n)
     &               - cellPhi(i-ii,j-jj,k-kk,l-ll,n)
            dphibarmin = min(abs(dphibarl),abs(dphibarr))
            if (dphifacemin.ge.dphibarmin) then
              dphichkl = dphifacel
              dphichkr = dphifacer
            else
              dphichkl = dphibarl
              dphichkr = dphibarr
            endif
            if (dphichkl*dphichkr.le.zero) then
               extremum = .true.
            else
               if (bigl) then
                  phil = -two * phir
               endif
               if (bigr) then
                  phir = -two * phil
               endif
            endif
        endif
        if (extremum) then
          ratio = D2PhiLim(i,j,k,l,n)
     &          / max(abs(D2Phi(i,j,k,l,n)),zeroTol)
          phil = phil * ratio
          phir = phir * ratio
        endif
        phiLeft(i+ii,j+jj,k+kk,l+ll,n)
     &    = phir + cellPhi(i,j,k,l,n);
        phiRight(i,j,k,l,n)
     &    = phil + cellPhi(i,j,k,l,n)
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine SELECTUPWIND_4D(
     &           facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,phiLeft
     &           ,iphiLeftlo0,iphiLeftlo1,iphiLeftlo2,iphiLeftlo3
     &           ,iphiLefthi0,iphiLefthi1,iphiLefthi2,iphiLefthi3
     &           ,nphiLeftcomp
     &           ,phiRight
     &           ,iphiRightlo0,iphiRightlo1,iphiRightlo2,iphiRightlo3
     &           ,iphiRighthi0,iphiRighthi1,iphiRighthi2,iphiRighthi3
     &           ,nphiRightcomp
     &           ,faceVel
     &           ,ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
     &           ,ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
     &           ,ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
     &           ,ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfacePhicomp
      integer ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
      integer ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
      REAL_T facePhi(
     &           ifacePhilo0:ifacePhihi0,
     &           ifacePhilo1:ifacePhihi1,
     &           ifacePhilo2:ifacePhihi2,
     &           ifacePhilo3:ifacePhihi3,
     &           0:nfacePhicomp-1)
      integer nphiLeftcomp
      integer iphiLeftlo0,iphiLeftlo1,iphiLeftlo2,iphiLeftlo3
      integer iphiLefthi0,iphiLefthi1,iphiLefthi2,iphiLefthi3
      REAL_T phiLeft(
     &           iphiLeftlo0:iphiLefthi0,
     &           iphiLeftlo1:iphiLefthi1,
     &           iphiLeftlo2:iphiLefthi2,
     &           iphiLeftlo3:iphiLefthi3,
     &           0:nphiLeftcomp-1)
      integer nphiRightcomp
      integer iphiRightlo0,iphiRightlo1,iphiRightlo2,iphiRightlo3
      integer iphiRighthi0,iphiRighthi1,iphiRighthi2,iphiRighthi3
      REAL_T phiRight(
     &           iphiRightlo0:iphiRighthi0,
     &           iphiRightlo1:iphiRighthi1,
     &           iphiRightlo2:iphiRighthi2,
     &           iphiRightlo3:iphiRighthi3,
     &           0:nphiRightcomp-1)
      integer ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
      integer ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
      REAL_T faceVel(
     &           ifaceVello0:ifaceVelhi0,
     &           ifaceVello1:ifaceVelhi1,
     &           ifaceVello2:ifaceVelhi2,
     &           ifaceVello3:ifaceVelhi3)
      integer ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
      integer ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
      integer i,j,k,l, n
      do n=0, nfacePhicomp-1
         
      do l = ifaceBoxlo3,ifaceBoxhi3
      do k = ifaceBoxlo2,ifaceBoxhi2
      do j = ifaceBoxlo1,ifaceBoxhi1
      do i = ifaceBoxlo0,ifaceBoxhi0

           if (faceVel(i,j,k,l).gt.zero) then
              facePhi(i,j,k,l,n)
     &          = phiLeft(i,j,k,l,n)
           else if (faceVel(i,j,k,l).lt.zero) then
              facePhi(i,j,k,l,n)
     &          = phiRight(i,j,k,l,n)
           else
              facePhi(i,j,k,l,n)
     &          = half * ( phiRight(i,j,k,l,n)
     &                   + phiLeft(i,j,k,l,n) )
           endif
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
