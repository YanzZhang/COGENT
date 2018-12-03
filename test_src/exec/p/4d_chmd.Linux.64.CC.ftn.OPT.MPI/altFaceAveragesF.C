#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
#define twentyseven (27.0d0)
#define fortyseven  (47.0d0)
#define eleven      (11.0d0)
#define thirteen    (13.0d0)
#define sixtieth    (1.0d0 / 60.0d0)
#define d0 (0.3d0)
#define d1 (0.6d0)
#define d2 (0.1d0)
#define thirteentwelfths (13.d0 / 12.d0)
#define quarter (0.25d0)
      subroutine BWENOFACEVALUES_4D(
     &           facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,faceVel
     &           ,ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
     &           ,ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
     &           ,ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
     &           ,ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
     &           ,idir
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
      integer ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
      integer ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
      REAL_T faceVel(
     &           ifaceVello0:ifaceVelhi0,
     &           ifaceVello1:ifaceVelhi1,
     &           ifaceVello2:ifaceVelhi2,
     &           ifaceVello3:ifaceVelhi3)
      integer ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
      integer ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
      integer idir
      integer n, i,j,k,l
      integer ii,jj,kk,ll
      REAL_T eps
      REAL_T fl,fr,bl,br,al,ar,wl,wr
      REAL_T c1l,c2l,c1r,c2r
      REAL_T wmax,wmin
      REAL_T onept5
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      onept5 = one + half
      eps = 1.d-6
      do n=0, (nfacePhicomp-1)
         
      do l = ifaceBoxlo3,ifaceBoxhi3
      do k = ifaceBoxlo2,ifaceBoxhi2
      do j = ifaceBoxlo1,ifaceBoxhi1
      do i = ifaceBoxlo0,ifaceBoxhi0

        fl = sixth*(
     &             - cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
     &             + five*cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &             + two*cellPhi(i,j,k,l,n) )
        fr = sixth*(
     &               two*cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &             + five*cellPhi(i,j,k,l,n)
     &             - cellPhi(i+ii,j+jj,k+kk,l+ll,n) )
        c1l = ( cellPhi(i,j,k,l,n)
     &        - two*cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &        + cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n) )
        c2l = ( cellPhi(i,j,k,l,n)
     &        - cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n) )
        c1r = ( cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &        - two*cellPhi(i,j,k,l,n)
     &        + cellPhi(i-ii,j-jj,k-kk,l-ll,n) )
        c2r = ( cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &        - cellPhi(i-ii,j-jj,k-kk,l-ll,n) )
        bl  = four*(c1l**2)*third+half*c1l*c2l+fourth*c2l**2
        br  = four*(c1r**2)*third-half*c1r*c2r+fourth*c2r**2
        al = one/((eps+bl)**2)
        ar = one/((eps+br)**2)
        wl = al/(al+ar)
        wr = ar/(al+ar)
        al = wl*(three*fourth+wl*(wl-onept5))
        ar = wr*(three*fourth+wr*(wr-onept5))
        wl = al/(al+ar)
        wr = ar/(al+ar)
        wmax = max(wl,wr)
        wmin = min(wl,wr)
        if( faceVel(i,j,k,l).gt.zero ) then
          wl = wmax
          wr = wmin
        else
          wl = wmin
          wr = wmax
        end if
        facePhi(i,j,k,l,n) = ( wl*fl + wr*fr )
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      function g( w, d )
      REAL_T g, w, d
      g = (d*(1-w)+(w-d)**2)/(d**2+w*(1-two*d))
      return
      end
      subroutine WENO5FACEVALUES_4D(
     &           facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,faceVel
     &           ,ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
     &           ,ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
     &           ,ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
     &           ,ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
     &           ,idir
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
      integer ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
      integer ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
      REAL_T faceVel(
     &           ifaceVello0:ifaceVelhi0,
     &           ifaceVello1:ifaceVelhi1,
     &           ifaceVello2:ifaceVelhi2,
     &           ifaceVello3:ifaceVelhi3)
      integer ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
      integer ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
      integer idir
      integer n, i,j,k,l
      integer ii,jj,kk,ll
      REAL_T eps
      REAL_T fl,fr,bl,br,al,ar,wl,wr
      REAL_T c0,c1,c2,c3
      REAL_T a0,a1,a2,asuminv
      REAL_T b0,b1,b2
      REAL_T v0,v1,v2
      REAL_T w0,w1,w2
      REAL_T g
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      eps = 1.d-6
      do n=0, (nfacePhicomp-1)
         
      do l = ifaceBoxlo3,ifaceBoxhi3
      do k = ifaceBoxlo2,ifaceBoxhi2
      do j = ifaceBoxlo1,ifaceBoxhi1
      do i = ifaceBoxlo0,ifaceBoxhi0

        if( faceVel(i,j,k,l).gt.zero ) then
          v0 = sixth*(
     &                 two*cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &               + five*cellPhi(i,j,k,l,n)
     &               - cellPhi(i+ii,j+jj,k+kk,l+ll,n) )
          v1 = sixth*(
     &               - cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
     &               + five*cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &               + two*cellPhi(i,j,k,l,n) )
          v2 = sixth*(
     &                 two*cellPhi(i-3*ii,j-3*jj,k-3*kk,l-3*ll,n)
     &               - seven*cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
     &               + eleven*cellPhi(i-ii,j-jj,k-kk,l-ll,n) )
          c0 = cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &       - cellPhi(i,j,k,l,n)
          c1 = cellPhi(i,j,k,l,n) 
     &       - cellPhi(i-ii,j-jj,k-kk,l-ll,n) 
          c2 = cellPhi(i-ii,j-jj,k-kk,l-ll,n) 
     &       - cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
          c3 = cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
     &       - cellPhi(i-3*ii,j-3*jj,k-3*kk,l-3*ll,n) 
        else
          v0 = sixth*(
     &                 two*cellPhi(i,j,k,l,n)
     &               + five*cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &               - cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n) )
          v1 = sixth*(
     &               - cellPhi(i+ii,j+jj,k+kk,l+ll,n) 
     &               + five*cellPhi(i,j,k,l,n)
     &               + two*cellPhi(i-ii,j-jj,k-kk,l-ll,n) )
          v2 = sixth*(
     &                 two*cellPhi(i+2*ii,j+2*jj,k+2*kk,l+2*ll,n)
     &               - seven*cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &               + eleven*cellPhi(i,j,k,l,n) )
          c0 = cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
     &       - cellPhi(i-ii,j-jj,k-kk,l-ll,n)
          c1 = cellPhi(i-ii,j-jj,k-kk,l-ll,n) 
     &       - cellPhi(i,j,k,l,n) 
          c2 = cellPhi(i,j,k,l,n) 
     &       - cellPhi(i+ii,j+jj,k+kk,l+ll,n)
          c3 = cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &       - cellPhi(i+2*ii,j+2*jj,k+2*kk,l+2*ll,n)
        end if
        b0 = thirteentwelfths*(c1-c0)**2+quarter*(three*c1-c0)**2
        b1 = thirteentwelfths*(c2-c1)**2+quarter*(c2+c1)**2
        b2 = thirteentwelfths*(c3-c2)**2+quarter*(c3-three*c2)**2
        a0 = d0/((eps+b0)**2)
        a1 = d1/((eps+b1)**2)
        a2 = d2/((eps+b2)**2)
        asuminv = one/(a0+a1+a2)
        w0 = a0*asuminv
        w1 = a1*asuminv
        w2 = a2*asuminv
        facePhi(i,j,k,l,n) = ( w0*v0 + w1*v1 + w2*v2 )
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine UW5FACEVALUES_4D(
     &           facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,faceVel
     &           ,ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
     &           ,ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
     &           ,ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
     &           ,ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
     &           ,idir
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
      integer ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
      integer ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
      REAL_T faceVel(
     &           ifaceVello0:ifaceVelhi0,
     &           ifaceVello1:ifaceVelhi1,
     &           ifaceVello2:ifaceVelhi2,
     &           ifaceVello3:ifaceVelhi3)
      integer ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
      integer ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
      integer idir
      integer n, i,j,k,l
      integer ii,jj,kk,ll
      REAL_T val
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      do n=0, (nfacePhicomp-1)
         
      do l = ifaceBoxlo3,ifaceBoxhi3
      do k = ifaceBoxlo2,ifaceBoxhi2
      do j = ifaceBoxlo1,ifaceBoxhi1
      do i = ifaceBoxlo0,ifaceBoxhi0

        if( faceVel(i,j,k,l).gt.zero ) then
          val = two         
     &          * cellPhi(i-3*ii,j-3*jj,k-3*kk,l-3*ll,n)
     &        - thirteen    
     &          * cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
     &        + fortyseven  
     &          * cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &        + twentyseven 
     &          * cellPhi(i,j,k,l,n)
     &        - three       
     &          * cellPhi(i+ii,j+jj,k+kk,l+ll,n)
        else
          val = two         
     &          * cellPhi(i+2*ii,j+2*jj,k+2*kk,l+2*ll,n)
     &        - thirteen    
     &          * cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &        + fortyseven  
     &          * cellPhi(i,j,k,l,n)
     &        + twentyseven 
     &          * cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &        - three       
     &          * cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
         end if
        facePhi(i,j,k,l,n) = val * sixtieth
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine UW3FACEVALUES_4D(
     &           facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,faceVel
     &           ,ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
     &           ,ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
     &           ,ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
     &           ,ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
     &           ,idir
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
      integer ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
      integer ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
      REAL_T faceVel(
     &           ifaceVello0:ifaceVelhi0,
     &           ifaceVello1:ifaceVelhi1,
     &           ifaceVello2:ifaceVelhi2,
     &           ifaceVello3:ifaceVelhi3)
      integer ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
      integer ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
      integer idir
      integer n, i,j,k,l
      integer ii,jj,kk,ll
      REAL_T val
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      do n=0, (nfacePhicomp-1)
         
      do l = ifaceBoxlo3,ifaceBoxhi3
      do k = ifaceBoxlo2,ifaceBoxhi2
      do j = ifaceBoxlo1,ifaceBoxhi1
      do i = ifaceBoxlo0,ifaceBoxhi0

        if( faceVel(i,j,k,l).gt.zero ) then
          val = - one         
     &          * cellPhi(i-2*ii,j-2*jj,k-2*kk,l-2*ll,n)
     &        + five    
     &          * cellPhi(i-ii,j-jj,k-kk,l-ll,n)
     &        + two
     &          * cellPhi(i,j,k,l,n)
        else
          val = - one         
     &          * cellPhi(i+ii,j+jj,k+kk,l+ll,n)
     &        + five  
     &          * cellPhi(i,j,k,l,n)
     &        + two
     &          * cellPhi(i-ii,j-jj,k-kk,l-ll,n)
         end if
        facePhi(i,j,k,l,n) = val * sixth
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine UW1FACEVALUES_4D(
     &           facePhi
     &           ,ifacePhilo0,ifacePhilo1,ifacePhilo2,ifacePhilo3
     &           ,ifacePhihi0,ifacePhihi1,ifacePhihi2,ifacePhihi3
     &           ,nfacePhicomp
     &           ,cellPhi
     &           ,icellPhilo0,icellPhilo1,icellPhilo2,icellPhilo3
     &           ,icellPhihi0,icellPhihi1,icellPhihi2,icellPhihi3
     &           ,ncellPhicomp
     &           ,faceVel
     &           ,ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
     &           ,ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
     &           ,ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
     &           ,ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
     &           ,idir
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
      integer ifaceVello0,ifaceVello1,ifaceVello2,ifaceVello3
      integer ifaceVelhi0,ifaceVelhi1,ifaceVelhi2,ifaceVelhi3
      REAL_T faceVel(
     &           ifaceVello0:ifaceVelhi0,
     &           ifaceVello1:ifaceVelhi1,
     &           ifaceVello2:ifaceVelhi2,
     &           ifaceVello3:ifaceVelhi3)
      integer ifaceBoxlo0,ifaceBoxlo1,ifaceBoxlo2,ifaceBoxlo3
      integer ifaceBoxhi0,ifaceBoxhi1,ifaceBoxhi2,ifaceBoxhi3
      integer idir
      integer n, i,j,k,l
      integer ii,jj,kk,ll
      REAL_T val
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      do n=0, (nfacePhicomp-1)
         
      do l = ifaceBoxlo3,ifaceBoxhi3
      do k = ifaceBoxlo2,ifaceBoxhi2
      do j = ifaceBoxlo1,ifaceBoxhi1
      do i = ifaceBoxlo0,ifaceBoxhi0

        if( faceVel(i,j,k,l).gt.zero ) then
          val = cellPhi(i-ii,j-jj,k-kk,l-ll,n)
        else
          val = cellPhi(i,j,k,l,n)
         end if
        facePhi(i,j,k,l,n) = val
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
