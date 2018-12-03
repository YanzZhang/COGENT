#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine FOURTH_ORDER_OUTFLOW_BC_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,vn
     &           ,ivnlo0,ivnlo1,ivnlo2
     &           ,ivnhi0,ivnhi1,ivnhi2
     &           ,nvncomp
     &           ,idir
     &           ,iside
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer nvncomp
      integer ivnlo0,ivnlo1,ivnlo2
      integer ivnhi0,ivnhi1,ivnhi2
      REAL_T vn(
     &           ivnlo0:ivnhi0,
     &           ivnlo1:ivnhi1,
     &           ivnlo2:ivnhi2,
     &           0:nvncomp-1)
      integer idir
      integer iside
      integer i0,i1,i2
      integer ii0,ii1,ii2
      integer iii0,iii1,iii2
      integer iv0,iv1,iv2
      integer ibeg0,ibeg1,ibeg2
      integer ng0,ng1,ng2
      integer isign
      integer n,d,ng(0:CH_SPACEDIM-1)
      isign = 2*iside-1
      
      iii0=isign*CHF_ID(0,idir)

      iii1=isign*CHF_ID(1,idir)

      iii2=isign*CHF_ID(2,idir)

      
      ng0 = ibdryboxhi0-ibdryboxlo0+1
      ng1 = ibdryboxhi1-ibdryboxlo1+1
      ng2 = ibdryboxhi2-ibdryboxlo2+1
      
      ibeg0 = ibdryboxlo0+CHF_ID(idir,0)*(1-iside)*ng0
      ibeg1 = ibdryboxlo1+CHF_ID(idir,1)*(1-iside)*ng1
      ibeg2 = ibdryboxlo2+CHF_ID(idir,2)*(1-iside)*ng2
      do n=0,nfcomp-1
      
      do i2 = ibdryboxlo2,ibdryboxhi2
      do i1 = ibdryboxlo1,ibdryboxhi1
      do i0 = ibdryboxlo0,ibdryboxhi0

        
        iv0 = i0+CHF_ID(idir,0)*(ibeg0-i0)
        iv1 = i1+CHF_ID(idir,1)*(ibeg1-i1)
        iv2 = i2+CHF_ID(idir,2)*(ibeg2-i2)
        if ( isign*vn(iv0,iv1,iv2,0).ge.zero ) then
          
          ii0 = i0+CHF_ID(idir,0)*(iside-1)*(2*(i0-ibeg0)+ng0+1)
          ii1 = i1+CHF_ID(idir,1)*(iside-1)*(2*(i1-ibeg1)+ng1+1)
          ii2 = i2+CHF_ID(idir,2)*(iside-1)*(2*(i2-ibeg2)+ng2+1)
          f(ii0,ii1,ii2,n) =
     &                   4*f(ii0-  iii0,ii1-  iii1,ii2-  iii2,n)
     &                 - 6*f(ii0-2*iii0,ii1-2*iii1,ii2-2*iii2,n)
     &                 + 4*f(ii0-3*iii0,ii1-3*iii1,ii2-3*iii2,n)
     &                 - 1*f(ii0-4*iii0,ii1-4*iii1,ii2-4*iii2,n)
        endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FOURTH_ORDER_NEUMANN_BC_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,ng
     &           ,idir
     &           ,iside
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer ng(0:2)
      integer idir
      integer iside
      integer i,j,k
      integer is,js,ks
      integer ibeg0,ibeg1,ibeg2
      integer n
      integer isign
      isign = 2*iside-1
      
      ibeg0 = ibdryboxlo0+CHF_ID(idir,0)*((1-iside)*(ng(0)+1)-1)
      ibeg1 = ibdryboxlo1+CHF_ID(idir,1)*((1-iside)*(ng(1)+1)-1)
      ibeg2 = ibdryboxlo2+CHF_ID(idir,2)*((1-iside)*(ng(2)+1)-1)
      do n=0,nfcomp-1
      
      do k = ibdryboxlo2,ibdryboxhi2
      do j = ibdryboxlo1,ibdryboxhi1
      do i = ibdryboxlo0,ibdryboxhi0

        
        is = i + CHF_ID(idir,0)*(2*(ibeg0-i)+isign)
        js = j + CHF_ID(idir,1)*(2*(ibeg1-j)+isign)
        ks = k + CHF_ID(idir,2)*(2*(ibeg2-k)+isign)
        f(i,j,k,n) = f(is,js,ks,n)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FOURTH_ORDER_DIRICHLET_BC_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,ng
     &           ,val
     &           ,ivallo0,ivallo1,ivallo2
     &           ,ivalhi0,ivalhi1,ivalhi2
     &           ,nvalcomp
     &           ,idir
     &           ,iside
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer ng(0:2)
      integer nvalcomp
      integer ivallo0,ivallo1,ivallo2
      integer ivalhi0,ivalhi1,ivalhi2
      REAL_T val(
     &           ivallo0:ivalhi0,
     &           ivallo1:ivalhi1,
     &           ivallo2:ivalhi2,
     &           0:nvalcomp-1)
      integer idir
      integer iside
      integer i0,i1,i2
      integer ii0,ii1,ii2
      integer iii0,iii1,iii2
      integer iv0,iv1,iv2
      integer ibeg0,ibeg1,ibeg2
      integer isign, gn(0:CH_SPACEDIM-1)
      integer n
      REAL_T thirteen
      parameter (thirteen=13.d0)
      isign = 2*iside-1
      
      iii0=isign*CHF_ID(0,idir)

      iii1=isign*CHF_ID(1,idir)

      iii2=isign*CHF_ID(2,idir)

      
      ibeg0 = ibdryboxlo0+CHF_ID(idir,0)*(1-iside)*ng(0)
      ibeg1 = ibdryboxlo1+CHF_ID(idir,1)*(1-iside)*ng(1)
      ibeg2 = ibdryboxlo2+CHF_ID(idir,2)*(1-iside)*ng(2)
      do n=0,nfcomp-1
      
      do i2 = ibdryboxlo2,ibdryboxhi2
      do i1 = ibdryboxlo1,ibdryboxhi1
      do i0 = ibdryboxlo0,ibdryboxhi0

        
        iv0 = i0+CHF_ID(idir,0)*(ibeg0-i0)
        iv1 = i1+CHF_ID(idir,1)*(ibeg1-i1)
        iv2 = i2+CHF_ID(idir,2)*(ibeg2-i2)
        
        ii0 = i0+CHF_ID(idir,0)*(iside-1)*(2*i0+ng(0)+1)
        ii1 = i1+CHF_ID(idir,1)*(iside-1)*(2*i1+ng(1)+1)
        ii2 = i2+CHF_ID(idir,2)*(iside-1)*(2*i2+ng(2)+1)
        
        gn(0) = ii0-ibeg0
        gn(1) = ii1-ibeg1
        gn(2) = ii2-ibeg2
        if (gn(idir).eq.1) then
          f(ii0,ii1,ii2,n) =
     *        third * (+ twelve   * val(iv0,iv1,iv2,n)
     &                 - thirteen * f(ii0-1*iii0,ii1-1*iii1,ii2-1*iii2,n)
     &                 + five     * f(ii0-2*iii0,ii1-2*iii1,ii2-2*iii2,n)
     &                 - one      * f(ii0-3*iii0,ii1-3*iii1,ii2-3*iii2,n))
        else if (gn(idir).eq.2) then
          f(ii0,ii1,ii2,n) =
     &                 + seven  * f(ii0-1*iii0,ii1-1*iii1,ii2-1*iii2,n)
     &                 - twelve * val(iv0,iv1,iv2,n)
     &                 + seven  * f(ii0-2*iii0,ii1-2*iii1,ii2-2*iii2,n)
     &                 - one    * f(ii0-3*iii0,ii1-3*iii1,ii2-3*iii2,n)
        else if (gn(idir).eq.3) then
          f(ii0,ii1,ii2,n) =
     &                 + five     * f(ii0-1*iii0,ii1-1*iii1,ii2-1*iii2,n)
     &                 - thirteen * f(ii0-2*iii0,ii1-2*iii1,ii2-2*iii2,n)
     *                 + twelve   * val(iv0,iv1,iv2,n)
     &                 - three    * f(ii0-3*iii0,ii1-3*iii1,ii2-3*iii2,n)
        else
          f(ii0,ii1,ii2,n) =
     &                   four * f(ii0-1*iii0,ii1-1*iii1,ii2-1*iii2,n)
     &                 - six  * f(ii0-2*iii0,ii1-2*iii1,ii2-2*iii2,n)
     &                 + four * f(ii0-3*iii0,ii1-3*iii1,ii2-3*iii2,n)
     &                 - one  * f(ii0-4*iii0,ii1-4*iii1,ii2-4*iii2,n)
        endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine SECOND_ORDER_DIRICHLET_BC_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iside
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer idir
      integer iside
      integer i,j,k
      integer ii,jj,kk
      integer i1,j1,k1
      integer iv0,iv1,iv2
      integer ibeg0,ibeg1,ibeg2
      integer ng0,ng1,ng2
      integer isign
      integer n,d,ng(0:CH_SPACEDIM-1),gn(0:CH_SPACEDIM-1)
      double precision val(0:nfcomp-1)
      isign = 2*iside-1
      
      i1 = isign * CHF_ID(idir,0)
      j1 = isign * CHF_ID(idir,1)
      k1 = isign * CHF_ID(idir,2)
      
      ng0 = ibdryboxhi0-ibdryboxlo0
      ng1 = ibdryboxhi1-ibdryboxlo1
      ng2 = ibdryboxhi2-ibdryboxlo2
      
      ibeg0 = ibdryboxlo0+CHF_ID(idir,0)*(1-iside)*ng0
      ibeg1 = ibdryboxlo1+CHF_ID(idir,1)*(1-iside)*ng1
      ibeg2 = ibdryboxlo2+CHF_ID(idir,2)*(1-iside)*ng2
      do n=0,nfcomp-1
      
      do k = ibdryboxlo2,ibdryboxhi2
      do j = ibdryboxlo1,ibdryboxhi1
      do i = ibdryboxlo0,ibdryboxhi0

        
        iv0 = i+CHF_ID(idir,0)*(ibeg0-i)
        iv1 = j+CHF_ID(idir,1)*(ibeg1-j)
        iv2 = k+CHF_ID(idir,2)*(ibeg2-k)
        val(n) = (three/two)*f(iv0,iv1,iv2,n) - (one/two)*f(iv0+i1,iv1+j1,iv2+k1,n)
        
        ii = i+CHF_ID(idir,0)*(iside-1)*(2*(i-ibeg0)+ng0)
        jj = j+CHF_ID(idir,1)*(iside-1)*(2*(j-ibeg1)+ng1)
        kk = k+CHF_ID(idir,2)*(iside-1)*(2*(k-ibeg2)+ng2)
        
        gn(0) = iabs(ii-ibeg0)
        gn(1) = iabs(jj-ibeg1)
        gn(2) = iabs(kk-ibeg2)
        f(ii,jj,kk,n) = two*(one + gn(idir)) * val(n) 
     &      - (one + two * gn(idir)) * f(iv0-i1,iv1-j1,iv2-k1,n)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
