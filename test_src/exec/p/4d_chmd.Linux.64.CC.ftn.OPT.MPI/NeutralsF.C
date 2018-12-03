#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine COMPUTE_IONIZATION_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2,irhslo3
     &           ,irhshi0,irhshi1,irhshi2,irhshi3
     &           ,fi
     &           ,ifilo0,ifilo1,ifilo2,ifilo3
     &           ,ifihi0,ifihi1,ifihi2,ifihi3
     &           ,nn
     &           ,innlo0,innlo1,innlo2,innlo3
     &           ,innhi0,innhi1,innhi2,innhi3
     &           ,sigmaV
     &           ,isigmaVlo0,isigmaVlo1,isigmaVlo2,isigmaVlo3
     &           ,isigmaVhi0,isigmaVhi1,isigmaVhi2,isigmaVhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer irhslo0,irhslo1,irhslo2,irhslo3
      integer irhshi0,irhshi1,irhshi2,irhshi3
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           irhslo3:irhshi3)
      integer ifilo0,ifilo1,ifilo2,ifilo3
      integer ifihi0,ifihi1,ifihi2,ifihi3
      REAL_T fi(
     &           ifilo0:ifihi0,
     &           ifilo1:ifihi1,
     &           ifilo2:ifihi2,
     &           ifilo3:ifihi3)
      integer innlo0,innlo1,innlo2,innlo3
      integer innhi0,innhi1,innhi2,innhi3
      REAL_T nn(
     &           innlo0:innhi0,
     &           innlo1:innhi1,
     &           innlo2:innhi2,
     &           innlo3:innhi3)
      integer isigmaVlo0,isigmaVlo1,isigmaVlo2,isigmaVlo3
      integer isigmaVhi0,isigmaVhi1,isigmaVhi2,isigmaVhi3
      REAL_T sigmaV(
     &           isigmaVlo0:isigmaVhi0,
     &           isigmaVlo1:isigmaVhi1,
     &           isigmaVlo2:isigmaVhi2,
     &           isigmaVlo3:isigmaVhi3)
      integer i,j,k,l
      
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==4
       rhs(i,j,k,l) = sigmaV(i,j,isigmaVlo2,isigmaVlo3)
     &                      * nn(i,j,innlo2,innlo3) 
     &                      * fi(i,j,k,l)
#else
       rhs(i,j,k,l) = sigmaV(i,j,k,isigmaVlo3,isigmaVlo4)
     &                      * nn(i,j,k,innlo3,innlo4) 
     &                      * fi(i,j,k,l)
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_CHARGE_EXCHANGE_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2,irhslo3
     &           ,irhshi0,irhshi1,irhshi2,irhshi3
     &           ,fi
     &           ,ifilo0,ifilo1,ifilo2,ifilo3
     &           ,ifihi0,ifihi1,ifihi2,ifihi3
     &           ,ni
     &           ,inilo0,inilo1,inilo2,inilo3
     &           ,inihi0,inihi1,inihi2,inihi3
     &           ,Ti
     &           ,iTilo0,iTilo1,iTilo2,iTilo3
     &           ,iTihi0,iTihi1,iTihi2,iTihi3
     &           ,nn
     &           ,innlo0,innlo1,innlo2,innlo3
     &           ,innhi0,innhi1,innhi2,innhi3
     &           ,fn
     &           ,ifnlo0,ifnlo1,ifnlo2,ifnlo3
     &           ,ifnhi0,ifnhi1,ifnhi2,ifnhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer irhslo0,irhslo1,irhslo2,irhslo3
      integer irhshi0,irhshi1,irhshi2,irhshi3
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           irhslo3:irhshi3)
      integer ifilo0,ifilo1,ifilo2,ifilo3
      integer ifihi0,ifihi1,ifihi2,ifihi3
      REAL_T fi(
     &           ifilo0:ifihi0,
     &           ifilo1:ifihi1,
     &           ifilo2:ifihi2,
     &           ifilo3:ifihi3)
      integer inilo0,inilo1,inilo2,inilo3
      integer inihi0,inihi1,inihi2,inihi3
      REAL_T ni(
     &           inilo0:inihi0,
     &           inilo1:inihi1,
     &           inilo2:inihi2,
     &           inilo3:inihi3)
      integer iTilo0,iTilo1,iTilo2,iTilo3
      integer iTihi0,iTihi1,iTihi2,iTihi3
      REAL_T Ti(
     &           iTilo0:iTihi0,
     &           iTilo1:iTihi1,
     &           iTilo2:iTihi2,
     &           iTilo3:iTihi3)
      integer innlo0,innlo1,innlo2,innlo3
      integer innhi0,innhi1,innhi2,innhi3
      REAL_T nn(
     &           innlo0:innhi0,
     &           innlo1:innhi1,
     &           innlo2:innhi2,
     &           innlo3:innhi3)
      integer ifnlo0,ifnlo1,ifnlo2,ifnlo3
      integer ifnhi0,ifnhi1,ifnhi2,ifnhi3
      REAL_T fn(
     &           ifnlo0:ifnhi0,
     &           ifnlo1:ifnhi1,
     &           ifnlo2:ifnhi2,
     &           ifnlo3:ifnhi3)
      integer i,j,k,l
      
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==4
       rhs(i,j,k,l) = ni(i,j,inilo2,inilo3) * fn(i,j,k,l)
     &                      - nn(i,j,inilo2,inilo3) * fi(i,j,k,l) 
       rhs(i,j,k,l) = rhs(i,j,k,l) * sqrt(Ti(i,j,iTilo2,iTilo3))
#else
       rhs(i,j,k,l) = ni(i,j,k,inilo3,inilo4) * fn(i,j,k,l)
     &                        - nn(i,j,k,inilo3,inilo4) * fi(i,j,k,l) 
       rhs(i,j,k,l) = rhs(i,j,k,l) * sqrt(Ti(i,j,k,iTilo3,iTilo4))
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_MODEL_CHARGE_EXCHANGE_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2,irhslo3
     &           ,irhshi0,irhshi1,irhshi2,irhshi3
     &           ,nn
     &           ,innlo0,innlo1,innlo2,innlo3
     &           ,innhi0,innhi1,innhi2,innhi3
     &           ,Ti
     &           ,iTilo0,iTilo1,iTilo2,iTilo3
     &           ,iTihi0,iTihi1,iTihi2,iTihi3
     &           ,u_par
     &           ,iu_parlo0,iu_parlo1,iu_parlo2,iu_parlo3
     &           ,iu_parhi0,iu_parhi1,iu_parhi2,iu_parhi3
     &           ,f0
     &           ,if0lo0,if0lo1,if0lo2,if0lo3
     &           ,if0hi0,if0hi1,if0hi2,if0hi3
     &           ,v_par
     &           ,iv_parlo0,iv_parlo1,iv_parlo2,iv_parlo3
     &           ,iv_parhi0,iv_parhi1,iv_parhi2,iv_parhi3
     &           ,mass
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer irhslo0,irhslo1,irhslo2,irhslo3
      integer irhshi0,irhshi1,irhshi2,irhshi3
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           irhslo3:irhshi3)
      integer innlo0,innlo1,innlo2,innlo3
      integer innhi0,innhi1,innhi2,innhi3
      REAL_T nn(
     &           innlo0:innhi0,
     &           innlo1:innhi1,
     &           innlo2:innhi2,
     &           innlo3:innhi3)
      integer iTilo0,iTilo1,iTilo2,iTilo3
      integer iTihi0,iTihi1,iTihi2,iTihi3
      REAL_T Ti(
     &           iTilo0:iTihi0,
     &           iTilo1:iTihi1,
     &           iTilo2:iTihi2,
     &           iTilo3:iTihi3)
      integer iu_parlo0,iu_parlo1,iu_parlo2,iu_parlo3
      integer iu_parhi0,iu_parhi1,iu_parhi2,iu_parhi3
      REAL_T u_par(
     &           iu_parlo0:iu_parhi0,
     &           iu_parlo1:iu_parhi1,
     &           iu_parlo2:iu_parhi2,
     &           iu_parlo3:iu_parhi3)
      integer if0lo0,if0lo1,if0lo2,if0lo3
      integer if0hi0,if0hi1,if0hi2,if0hi3
      REAL_T f0(
     &           if0lo0:if0hi0,
     &           if0lo1:if0hi1,
     &           if0lo2:if0hi2,
     &           if0lo3:if0hi3)
      integer iv_parlo0,iv_parlo1,iv_parlo2,iv_parlo3
      integer iv_parhi0,iv_parhi1,iv_parhi2,iv_parhi3
      REAL_T v_par(
     &           iv_parlo0:iv_parhi0,
     &           iv_parlo1:iv_parhi1,
     &           iv_parlo2:iv_parhi2,
     &           iv_parlo3:iv_parhi3)
      REAL_T mass
      integer i,j,k,l
      
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==4
       rhs(i,j,k,l) = u_par(i,j,iu_parlo2,iu_parlo3) * f0(i,j,k,l)
     &                      * v_par(i,j,k,l) * mass / Ti(i,j,iTilo2,iTilo3)
       rhs(i,j,k,l) = rhs(i,j,k,l) * sqrt(Ti(i,j,iTilo2,iTilo3))
     &                      * nn(i,j,innlo2,innlo3)
#else
       rhs(i,j,k,l) = u_par(i,j,k,iu_parlo3,iu_parlo4) * f0(i,j,k,l)
     &                        * v_par(i,j,k,l) * mass / Ti(i,j,k,iTilo3,iTilo4)
       rhs(i,j,k,l) = rhs(i,j,k,l) * sqrt(Ti(i,j,k,iTilo3,iTilo4))
     &                        * nn(i,j,k,innlo3,innlo4)
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
