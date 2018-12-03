#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine INCREMENTFACEPROD_4D(
     &           prod
     &           ,iprodlo0,iprodlo1,iprodlo2,iprodlo3
     &           ,iprodhi0,iprodhi1,iprodhi2,iprodhi3
     &           ,nprodcomp
     &           ,u
     &           ,iulo0,iulo1,iulo2,iulo3
     &           ,iuhi0,iuhi1,iuhi2,iuhi3
     &           ,nucomp
     &           ,v
     &           ,ivlo0,ivlo1,ivlo2,ivlo3
     &           ,ivhi0,ivhi1,ivhi2,ivhi3
     &           ,nvcomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nprodcomp
      integer iprodlo0,iprodlo1,iprodlo2,iprodlo3
      integer iprodhi0,iprodhi1,iprodhi2,iprodhi3
      REAL_T prod(
     &           iprodlo0:iprodhi0,
     &           iprodlo1:iprodhi1,
     &           iprodlo2:iprodhi2,
     &           iprodlo3:iprodhi3,
     &           0:nprodcomp-1)
      integer nucomp
      integer iulo0,iulo1,iulo2,iulo3
      integer iuhi0,iuhi1,iuhi2,iuhi3
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2,
     &           iulo3:iuhi3,
     &           0:nucomp-1)
      integer nvcomp
      integer ivlo0,ivlo1,ivlo2,ivlo3
      integer ivhi0,ivhi1,ivhi2,ivhi3
      REAL_T v(
     &           ivlo0:ivhi0,
     &           ivlo1:ivhi1,
     &           ivlo2:ivhi2,
     &           ivlo3:ivhi3,
     &           0:nvcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer i,j,k,l, n
      integer d, nn
      do n=0, nucomp-1
        do d=0, CH_SPACEDIM-1
          nn = n + d * nucomp
          
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          prod(i,j,k,l,nn) = prod(i,j,k,l,nn)
     &             +u(i,j,k,l,n)*v(i,j,k,l,d)
          
      enddo
      enddo
      enddo
      enddo
       enddo
      enddo
      return
      end
      subroutine INCREMENTFACEPRODNORMAL_4D(
     &           prod
     &           ,iprodlo0,iprodlo1,iprodlo2,iprodlo3
     &           ,iprodhi0,iprodhi1,iprodhi2,iprodhi3
     &           ,nprodcomp
     &           ,u
     &           ,iulo0,iulo1,iulo2,iulo3
     &           ,iuhi0,iuhi1,iuhi2,iuhi3
     &           ,nucomp
     &           ,v
     &           ,ivlo0,ivlo1,ivlo2,ivlo3
     &           ,ivhi0,ivhi1,ivhi2,ivhi3
     &           ,nvcomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nprodcomp
      integer iprodlo0,iprodlo1,iprodlo2,iprodlo3
      integer iprodhi0,iprodhi1,iprodhi2,iprodhi3
      REAL_T prod(
     &           iprodlo0:iprodhi0,
     &           iprodlo1:iprodhi1,
     &           iprodlo2:iprodhi2,
     &           iprodlo3:iprodhi3,
     &           0:nprodcomp-1)
      integer nucomp
      integer iulo0,iulo1,iulo2,iulo3
      integer iuhi0,iuhi1,iuhi2,iuhi3
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2,
     &           iulo3:iuhi3,
     &           0:nucomp-1)
      integer nvcomp
      integer ivlo0,ivlo1,ivlo2,ivlo3
      integer ivhi0,ivhi1,ivhi2,ivhi3
      REAL_T v(
     &           ivlo0:ivhi0,
     &           ivlo1:ivhi1,
     &           ivlo2:ivhi2,
     &           ivlo3:ivhi3,
     &           0:nvcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer i,j,k,l, n
      do n=0, nucomp-1
         
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

         prod(i,j,k,l,n) = prod(i,j,k,l,n)
     &        +u(i,j,k,l,n)*v(i,j,k,l,n)
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine INCREMENTFACEPRODGRAD_4D(
     &           prod
     &           ,iprodlo0,iprodlo1,iprodlo2,iprodlo3
     &           ,iprodhi0,iprodhi1,iprodhi2,iprodhi3
     &           ,nprodcomp
     &           ,u
     &           ,iulo0,iulo1,iulo2,iulo3
     &           ,iuhi0,iuhi1,iuhi2,iuhi3
     &           ,nucomp
     &           ,v
     &           ,ivlo0,ivlo1,ivlo2,ivlo3
     &           ,ivhi0,ivhi1,ivhi2,ivhi3
     &           ,nvcomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,dx
     &           ,factor
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nprodcomp
      integer iprodlo0,iprodlo1,iprodlo2,iprodlo3
      integer iprodhi0,iprodhi1,iprodhi2,iprodhi3
      REAL_T prod(
     &           iprodlo0:iprodhi0,
     &           iprodlo1:iprodhi1,
     &           iprodlo2:iprodhi2,
     &           iprodlo3:iprodhi3,
     &           0:nprodcomp-1)
      integer nucomp
      integer iulo0,iulo1,iulo2,iulo3
      integer iuhi0,iuhi1,iuhi2,iuhi3
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2,
     &           iulo3:iuhi3,
     &           0:nucomp-1)
      integer nvcomp
      integer ivlo0,ivlo1,ivlo2,ivlo3
      integer ivhi0,ivhi1,ivhi2,ivhi3
      REAL_T v(
     &           ivlo0:ivhi0,
     &           ivlo1:ivhi1,
     &           ivlo2:ivhi2,
     &           ivlo3:ivhi3,
     &           0:nvcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      REAL_T dx
      REAL_T factor
      integer dir
      integer i,j,k,l, n, nn, d
      integer ii,jj,kk,ll
      REAL_T qtrOnDx2, du, dv
      ii = CHF_ID(0,dir)
                jj = CHF_ID(1,dir)
                kk = CHF_ID(2,dir)
                ll = CHF_ID(3,dir)
      qtrOnDx2 = (half/dx)**2;
      do n=0, nucomp-1
        do d=0, CH_SPACEDIM-1
          nn = n + d * nucomp
          
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            du = (u(i+ii,j+jj,k+kk,l+ll,n) - u(i-ii,j-jj,k-kk,l-ll,n))
            dv = (v(i+ii,j+jj,k+kk,l+ll,d) - v(i-ii,j-jj,k-kk,l-ll,d))
            prod(i,j,k,l,nn) = prod(i,j,k,l,nn)
     &                             +factor*qtrOnDx2*du*dv
         
      enddo
      enddo
      enddo
      enddo
        enddo
      enddo
      return
      end
      subroutine INCREMENTFACEPRODGRADNORMAL_4D(
     &           prod
     &           ,iprodlo0,iprodlo1,iprodlo2,iprodlo3
     &           ,iprodhi0,iprodhi1,iprodhi2,iprodhi3
     &           ,nprodcomp
     &           ,u
     &           ,iulo0,iulo1,iulo2,iulo3
     &           ,iuhi0,iuhi1,iuhi2,iuhi3
     &           ,nucomp
     &           ,v
     &           ,ivlo0,ivlo1,ivlo2,ivlo3
     &           ,ivhi0,ivhi1,ivhi2,ivhi3
     &           ,nvcomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,dx
     &           ,factor
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nprodcomp
      integer iprodlo0,iprodlo1,iprodlo2,iprodlo3
      integer iprodhi0,iprodhi1,iprodhi2,iprodhi3
      REAL_T prod(
     &           iprodlo0:iprodhi0,
     &           iprodlo1:iprodhi1,
     &           iprodlo2:iprodhi2,
     &           iprodlo3:iprodhi3,
     &           0:nprodcomp-1)
      integer nucomp
      integer iulo0,iulo1,iulo2,iulo3
      integer iuhi0,iuhi1,iuhi2,iuhi3
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2,
     &           iulo3:iuhi3,
     &           0:nucomp-1)
      integer nvcomp
      integer ivlo0,ivlo1,ivlo2,ivlo3
      integer ivhi0,ivhi1,ivhi2,ivhi3
      REAL_T v(
     &           ivlo0:ivhi0,
     &           ivlo1:ivhi1,
     &           ivlo2:ivhi2,
     &           ivlo3:ivhi3,
     &           0:nvcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      REAL_T dx
      REAL_T factor
      integer dir
      integer i,j,k,l, n
      integer ii,jj,kk,ll
      REAL_T qtrOnDx2, du, dv
      ii = CHF_ID(0,dir)
                jj = CHF_ID(1,dir)
                kk = CHF_ID(2,dir)
                ll = CHF_ID(3,dir)
      qtrOnDx2 = (half/dx)**2;
      do n=0, nucomp-1
         
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            du = (u(i+ii,j+jj,k+kk,l+ll,n) - u(i-ii,j-jj,k-kk,l-ll,n))
            dv = (v(i+ii,j+jj,k+kk,l+ll,0) - v(i-ii,j-jj,k-kk,l-ll,0))
            prod(i,j,k,l,n) = prod(i,j,k,l,n)
     &                             +factor*qtrOnDx2*du*dv
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
