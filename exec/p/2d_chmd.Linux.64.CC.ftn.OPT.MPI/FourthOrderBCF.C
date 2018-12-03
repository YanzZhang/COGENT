#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine FOURTH_ORDER_OUTFLOW_BC_2D(
     &           f
     &           ,iflo0,iflo1
     &           ,ifhi0,ifhi1
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1
     &           ,ibdryboxhi0,ibdryboxhi1
     &           ,vn
     &           ,ivnlo0,ivnlo1
     &           ,ivnhi0,ivnhi1
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
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1
      integer ibdryboxhi0,ibdryboxhi1
      integer nvncomp
      integer ivnlo0,ivnlo1
      integer ivnhi0,ivnhi1
      REAL_T vn(
     &           ivnlo0:ivnhi0,
     &           ivnlo1:ivnhi1,
     &           0:nvncomp-1)
      integer idir
      integer iside
      integer i0,i1
      integer ii0,ii1
      integer iii0,iii1
      integer iv0,iv1
      integer ibeg0,ibeg1
      integer ng0,ng1
      integer isign
      integer n,d,ng(0:CH_SPACEDIM-1)
      isign = 2*iside-1
      
      iii0=isign*CHF_ID(0,idir)

      iii1=isign*CHF_ID(1,idir)

      
      ng0 = ibdryboxhi0-ibdryboxlo0+1
      ng1 = ibdryboxhi1-ibdryboxlo1+1
      
      ibeg0 = ibdryboxlo0+CHF_ID(idir,0)*(1-iside)*ng0
      ibeg1 = ibdryboxlo1+CHF_ID(idir,1)*(1-iside)*ng1
      do n=0,nfcomp-1
      
      do i1 = ibdryboxlo1,ibdryboxhi1
      do i0 = ibdryboxlo0,ibdryboxhi0

        
        iv0 = i0+CHF_ID(idir,0)*(ibeg0-i0)
        iv1 = i1+CHF_ID(idir,1)*(ibeg1-i1)
        if ( isign*vn(iv0,iv1,0).ge.zero ) then
          
          ii0 = i0+CHF_ID(idir,0)*(iside-1)*(2*(i0-ibeg0)+ng0+1)
          ii1 = i1+CHF_ID(idir,1)*(iside-1)*(2*(i1-ibeg1)+ng1+1)
          f(ii0,ii1,n) =
     &                   4*f(ii0-  iii0,ii1-  iii1,n)
     &                 - 6*f(ii0-2*iii0,ii1-2*iii1,n)
     &                 + 4*f(ii0-3*iii0,ii1-3*iii1,n)
     &                 - 1*f(ii0-4*iii0,ii1-4*iii1,n)
        endif
      
      enddo
      enddo
      enddo
      return
      end
      subroutine FOURTH_ORDER_NEUMANN_BC_2D(
     &           f
     &           ,iflo0,iflo1
     &           ,ifhi0,ifhi1
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1
     &           ,ibdryboxhi0,ibdryboxhi1
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
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1
      integer ibdryboxhi0,ibdryboxhi1
      integer ng(0:1)
      integer idir
      integer iside
      integer i,j
      integer is,js
      integer ibeg0,ibeg1
      integer n
      integer isign
      isign = 2*iside-1
      
      ibeg0 = ibdryboxlo0+CHF_ID(idir,0)*((1-iside)*(ng(0)+1)-1)
      ibeg1 = ibdryboxlo1+CHF_ID(idir,1)*((1-iside)*(ng(1)+1)-1)
      do n=0,nfcomp-1
      
      do j = ibdryboxlo1,ibdryboxhi1
      do i = ibdryboxlo0,ibdryboxhi0

        
        is = i + CHF_ID(idir,0)*(2*(ibeg0-i)+isign)
        js = j + CHF_ID(idir,1)*(2*(ibeg1-j)+isign)
        f(i,j,n) = f(is,js,n)
      
      enddo
      enddo
      enddo
      return
      end
      subroutine FOURTH_ORDER_DIRICHLET_BC_2D(
     &           f
     &           ,iflo0,iflo1
     &           ,ifhi0,ifhi1
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1
     &           ,ibdryboxhi0,ibdryboxhi1
     &           ,ng
     &           ,val
     &           ,ivallo0,ivallo1
     &           ,ivalhi0,ivalhi1
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
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1
      integer ibdryboxhi0,ibdryboxhi1
      integer ng(0:1)
      integer nvalcomp
      integer ivallo0,ivallo1
      integer ivalhi0,ivalhi1
      REAL_T val(
     &           ivallo0:ivalhi0,
     &           ivallo1:ivalhi1,
     &           0:nvalcomp-1)
      integer idir
      integer iside
      integer i0,i1
      integer ii0,ii1
      integer iii0,iii1
      integer iv0,iv1
      integer ibeg0,ibeg1
      integer isign, gn(0:CH_SPACEDIM-1)
      integer n
      REAL_T thirteen
      parameter (thirteen=13.d0)
      isign = 2*iside-1
      
      iii0=isign*CHF_ID(0,idir)

      iii1=isign*CHF_ID(1,idir)

      
      ibeg0 = ibdryboxlo0+CHF_ID(idir,0)*(1-iside)*ng(0)
      ibeg1 = ibdryboxlo1+CHF_ID(idir,1)*(1-iside)*ng(1)
      do n=0,nfcomp-1
      
      do i1 = ibdryboxlo1,ibdryboxhi1
      do i0 = ibdryboxlo0,ibdryboxhi0

        
        iv0 = i0+CHF_ID(idir,0)*(ibeg0-i0)
        iv1 = i1+CHF_ID(idir,1)*(ibeg1-i1)
        
        ii0 = i0+CHF_ID(idir,0)*(iside-1)*(2*i0+ng(0)+1)
        ii1 = i1+CHF_ID(idir,1)*(iside-1)*(2*i1+ng(1)+1)
        
        gn(0) = ii0-ibeg0
        gn(1) = ii1-ibeg1
        if (gn(idir).eq.1) then
          f(ii0,ii1,n) =
     *        third * (+ twelve   * val(iv0,iv1,n)
     &                 - thirteen * f(ii0-1*iii0,ii1-1*iii1,n)
     &                 + five     * f(ii0-2*iii0,ii1-2*iii1,n)
     &                 - one      * f(ii0-3*iii0,ii1-3*iii1,n))
        else if (gn(idir).eq.2) then
          f(ii0,ii1,n) =
     &                 + seven  * f(ii0-1*iii0,ii1-1*iii1,n)
     &                 - twelve * val(iv0,iv1,n)
     &                 + seven  * f(ii0-2*iii0,ii1-2*iii1,n)
     &                 - one    * f(ii0-3*iii0,ii1-3*iii1,n)
        else if (gn(idir).eq.3) then
          f(ii0,ii1,n) =
     &                 + five     * f(ii0-1*iii0,ii1-1*iii1,n)
     &                 - thirteen * f(ii0-2*iii0,ii1-2*iii1,n)
     *                 + twelve   * val(iv0,iv1,n)
     &                 - three    * f(ii0-3*iii0,ii1-3*iii1,n)
        else
          f(ii0,ii1,n) =
     &                   four * f(ii0-1*iii0,ii1-1*iii1,n)
     &                 - six  * f(ii0-2*iii0,ii1-2*iii1,n)
     &                 + four * f(ii0-3*iii0,ii1-3*iii1,n)
     &                 - one  * f(ii0-4*iii0,ii1-4*iii1,n)
        endif
      
      enddo
      enddo
      enddo
      return
      end
      subroutine SECOND_ORDER_DIRICHLET_BC_2D(
     &           f
     &           ,iflo0,iflo1
     &           ,ifhi0,ifhi1
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1
     &           ,ibdryboxhi0,ibdryboxhi1
     &           ,idir
     &           ,iside
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfcomp
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1
      integer ibdryboxhi0,ibdryboxhi1
      integer idir
      integer iside
      integer i,j
      integer ii,jj
      integer i1,j1
      integer iv0,iv1
      integer ibeg0,ibeg1
      integer ng0,ng1
      integer isign
      integer n,d,ng(0:CH_SPACEDIM-1),gn(0:CH_SPACEDIM-1)
      double precision val(0:nfcomp-1)
      isign = 2*iside-1
      
      i1 = isign * CHF_ID(idir,0)
      j1 = isign * CHF_ID(idir,1)
      
      ng0 = ibdryboxhi0-ibdryboxlo0
      ng1 = ibdryboxhi1-ibdryboxlo1
      
      ibeg0 = ibdryboxlo0+CHF_ID(idir,0)*(1-iside)*ng0
      ibeg1 = ibdryboxlo1+CHF_ID(idir,1)*(1-iside)*ng1
      do n=0,nfcomp-1
      
      do j = ibdryboxlo1,ibdryboxhi1
      do i = ibdryboxlo0,ibdryboxhi0

        
        iv0 = i+CHF_ID(idir,0)*(ibeg0-i)
        iv1 = j+CHF_ID(idir,1)*(ibeg1-j)
        val(n) = (three/two)*f(iv0,iv1,n) - (one/two)*f(iv0+i1,iv1+j1,n)
        
        ii = i+CHF_ID(idir,0)*(iside-1)*(2*(i-ibeg0)+ng0)
        jj = j+CHF_ID(idir,1)*(iside-1)*(2*(j-ibeg1)+ng1)
        
        gn(0) = iabs(ii-ibeg0)
        gn(1) = iabs(jj-ibeg1)
        f(ii,jj,n) = two*(one + gn(idir)) * val(n) 
     &      - (one + two * gn(idir)) * f(iv0-i1,iv1-j1,n)
      
      enddo
      enddo
      enddo
      return
      end
